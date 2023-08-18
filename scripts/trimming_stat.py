import pandas as pd
import sys
import os
import re
def read_fastqTr_metric(infile):
    ''' parse trimming statistics, *'.umi.r1.fqTr.metrics' '''
    with open(infile, 'r') as f:
        start_read = False
        summary_lines = []
        for line in f:
            if '=== Summary ===' in line:
                start_read = True
            if 'Adapter' in line:
                start_read = False
                break
            
            if start_read and len(line.rstrip())>0 and 'Summary' not in line and '==' not in line:
                line = line.rstrip()
                
                item, values = re.split(r':\s+', line)
                item = item.replace(':', '')
                
                if '%' in values:
                    number, perc = values.split('(')
                else:
                    number = values
                    perc = None
                
                number = int(number.replace(',', '').replace(' bp', ''))
                if perc is not None:
                    perc = float(perc.replace('%', '').replace(')', ''))
                    
                
                summary_lines.append([item, number, perc])
                
    # sometimes duplicated names
    df = pd.DataFrame(summary_lines, columns = ['item', 'number', 'percentage'])
    
    for index, row in df.iterrows():
        if 'Total basepairs processed' in row['item']:
            col = 'Total basepairs processed'
        if 'Total written (filtered)' in row['item']:
            col = 'Total written (filtered)'
        if 'Quality-trimmed' in row['item']:
            col = 'Quality trimmed'
        if 'Read 1' in row['item'] or 'Read 2' in row['item'] and 'adapter' not in row['item']:
            try:
                df.loc[index, 'item'] = col + row['item'] # rename so that it does not duplicate
            except:
                print(row['item'])
                # change name
    
    
    return df.set_index('item')

if __name__=='__main__':
    file_list = sys.argv[1].split(' ')
    print(file_list)
    out = sys.argv[2]

    ndf = []
    pdf = []
    names = []
    for file in file_list:
        name = os.path.basename(file)
        s = read_fastqTr_metric(file)
        ndf.append(s['number'])
        pdf.append(s['percentage'])
        names.append(name)
    print(sorted(names))
    print(len(names), len(set(names)))
    print('INDEX:', s.index)
    print('INDEX:', sorted(s.index.tolist()))
    ndf = pd.DataFrame(ndf, index = names)
    pdf = pd.DataFrame(pdf, index = names)

    pdf.dropna(axis = 'columns', how = 'all', inplace = True)
    pdf.columns = ['% ' + c for c in pdf.columns]

    total = pd.concat([ndf, pdf], axis = 1)

    total.to_csv(out)
