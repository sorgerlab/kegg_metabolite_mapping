import csv
import re
import sys
import lxml.etree
import pandas as pd


def main(argv=sys.argv):

    global hmdb, metabolites, matches

    args = argv[1:]
    if len(args) != 1:
        print("Usage: map_metabolites.py input.csv > output.csv")
        return
    in_path = args[0]

    ## Parse metabolite names from input file column headers.
    with open(in_path, 'r') as f:
        reader = csv.reader(f)
        headers = next(reader)
    # The input format has three columns of row headers. We need to skip over
    # those to get to the column headers.
    metabolites = headers[3:]
    # Each header has the word "Results" tacked on the end, so strip that.
    metabolites = [re.sub(r' Results$', '', m) for m in metabolites]
    metabolites = pd.DataFrame({'name': metabolites})

    ## Parse HMDB metabolite XML file.
    context = lxml.etree.iterparse('input/hmdb_metabolites.xml',
                                   tag='{http://www.hmdb.ca}metabolite')
    ns = {'m': 'http://www.hmdb.ca'}
    # Build parallel lists of names and database IDs, then turn them into a
    # DataFrame.
    names = []
    kegg_ids = []
    hmdb_ids = []
    for i, (action, m_elt) in enumerate(context):
        # Extract and clean up values.
        kegg_id = m_elt.findtext('m:kegg_id', '', namespaces=ns).strip()
        if kegg_id == '':
            kegg_id = None
        hmdb_id = m_elt.findtext('m:accession', namespaces=ns).strip()
        name = m_elt.findtext('m:name', namespaces=ns).strip()
        synonym_elts = m_elt.findall('m:synonyms/m:synonym', namespaces=ns)
        # Append values to lists.
        n_names = 1 + len(synonym_elts)
        hmdb_ids += [hmdb_id] * n_names
        kegg_ids += [kegg_id] * n_names
        names += [name] + [e.text for e in synonym_elts]
        m_elt.clear()
    hmdb = pd.DataFrame(
        {'hmdb_id': hmdb_ids, 'kegg_id': kegg_ids, 'name': names}
    )
    # This DataFrame ends up with some redundant rows as HMDB sometimes includes
    # the primary name in the synonyms list.
    hmdb.drop_duplicates(inplace=True)

    # Delete some overeager synonyms that lead to collisions in our data.
    for hid, name in [('HMDB0003192', "AICAR"),
                      ('HMDB0012305', "Uridine 5'-diphosphogalactose")]:
        idx = (hmdb['hmdb_id'] == hid) & (hmdb['name'] == name)
        hmdb.loc[idx, 'name'] = None

    # Generate an extra normalized name column to aid in string matching.
    for df in metabolites, hmdb:
        df['name_normalized'] = (df['name'].str.replace(r"[ ',_-]", "")
                                 .str.lower())
    assert metabolites.name_normalized.is_unique, "Name normalization collision"
    # Filter out records with duplicate normalized names for the same kegg_id.
    # Leave records with no kegg_id alone.
    hmdb = pd.concat([
        hmdb[hmdb.kegg_id.notnull()].drop_duplicates(['kegg_id',
                                                      'name_normalized']),
        hmdb[hmdb.kegg_id.isnull()].drop_duplicates(['hmdb_id',
                                                    'name_normalized'])
    ])

    matches = pd.merge(metabolites, hmdb, how='left', suffixes=['', '_hmdb'],
                     on='name_normalized')

    # Append manually curated entries.
    matches['manual'] = False
    manual = pd.read_csv('resources/curated_metabolites.csv')
    manual['manual'] = True
    matches = matches.append(manual, ignore_index=True)

    # If duplicate matches for a given metabolite name, throw away any
    # non-manual entries that have no kegg_id.
    discard_idx = (
        matches.duplicated(['name'], keep=False) & matches['kegg_id'].isnull() &
        ~matches['manual']
    )
    matches = matches[~discard_idx].reset_index(drop=True)

    # Reorder matches by original input file order.
    matches = matches.set_index('name').loc[metabolites['name']].reset_index()

    # Ensure all metabolites are represented.
    assert len(matches['name']) == len(metabolites['name'])
    assert set(matches['name']) == set(metabolites['name'])

    # Display mapping table.
    out = matches[['name','kegg_id','hmdb_id']].fillna('')
    out.to_csv(sys.stdout, index=False)


if __name__ == '__main__':
    main()
