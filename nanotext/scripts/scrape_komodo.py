import json

from bs4 import BeautifulSoup
import requests
from tqdm import tqdm


def parse_html_table(table, includes_header=True):
    '''
    https://srome.github.io/Parsing-HTML-Tables-in-Python-with-BeautifulSoup-and-pandas/
    
    Example:

    url = 'http://komodo.modelseed.org/servlet/KomodoTomcatServerSideUtilitiesModelSeed?MediaInfo='
    num = 3094
    r = requests.get(f'{url}{num}')
    soup = BeautifulSoup(r.text)
    tables = soup.find_all('table')
    # table 2 always seems to be what we are after
    parse_html_table(tables[2])
    '''
    import pandas as pd

    n_columns = 0
    n_rows=0
    column_names = []
    
    # Find number of rows and columns
    # we also find the column titles if we can
    for row in table.find_all('tr'):
        # Determine the number of rows in the table
        td_tags = row.find_all('td')
        
        if len(td_tags) > 0:
            n_rows += 1

            if n_columns == 0:
                # Set the number of columns for our table
                n_columns = len(td_tags)

        # Handle column names if we find them
        th_tags = row.find_all('th') 
        if len(th_tags) > 0 and len(column_names) == 0:
            for th in th_tags:
                column_names.append(th.get_text())

    # Safeguard on Column Titles
    if len(column_names) > 0 and len(column_names) != n_columns:
        raise Exception("Column titles do not match the number of columns")
    
    columns = column_names if len(column_names) > 0 else range(0, n_columns)
    df = pd.DataFrame(columns=columns, index=range(0, n_rows))
    
    row_marker = 0
    for row in table.find_all('tr'):
        column_marker = 0
        columns = row.find_all('td')
        
        for column in columns:
            df.iat[row_marker,column_marker] = column.get_text()
            column_marker += 1
        
        if len(columns) > 0:
            row_marker += 1
    
    # Convert to float if possible
    for col in df:
        try:
            df[col] = df[col].astype(float)
        except ValueError:
            pass
    
    if not column_names and includes_header:
        # stackoverflow.com/questions/26147180
        names = df.iloc[0]
        # print('No table header marked in html, will use the first row ...')
        df = pd.DataFrame(df.values[1:], columns=names)

    return df


# get list of all media
url = 'http://komodo.modelseed.org/servlet/KomodoTomcatServerSideUtilitiesModelSeed?MediaList'
r = requests.get(url)
soup = BeautifulSoup(r.text)
# the media table is the 2nd table on the webpage
df_media = parse_html_table(soup.find_all('table')[1])
# df_media.head().iloc[0]['Instructions']
# 'http://www.dsmz.de/microorganisms/medium/pdf/DSMZ_Medium1.pdf'
# We could actually use this as UID of the DSMZ recipe or as document tag.
doctags = []
for _, i in df_media.iterrows():
    doctags.append(
        i['Instructions'].strip('.pdf').split('/')[-1].strip('DSMZ_Medium'))
df_media['doctags'] = doctags

# save
fp = '/Users/phi/data_local/bacdive_media/media.tsv'
df_media.to_csv(fp, sep='\t', index=False)


# parse individual media
url = 'http://komodo.modelseed.org/servlet/KomodoTomcatServerSideUtilitiesModelSeed?MediaInfo='
fp = '/path/to/bacdive_media/media/'

for _, i in tqdm(df_media.iterrows()):
    r = requests.get(f'{url}{i["ID"]}')
    soup = BeautifulSoup(r.text)
    medium = parse_html_table(soup.find_all('table')[2])
    medium.to_csv(f'{fp}{i["ID"]}.tsv', sep='\t', index=False)  # save


# get medium-strain mapping
url = 'http://komodo.modelseed.org/servlet/KomodoTomcatServerSideUtilitiesModelSeed?OrganismMedia'
r = requests.get(url)
soup = BeautifulSoup(r.text)
df_bugs = parse_html_table(soup.find_all('table')[1])

media = []
for _, i in df_bugs.iterrows():
    media.append(i['Media list'].strip(' \t\t\t\t'))
df_bugs['Media list'] = media
fp = '/path/to/bacdive_media/bugs.tsv'
df_bugs.to_csv(fp, sep='\t', index=False)
