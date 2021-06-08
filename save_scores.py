from openpyxl import load_workbook
import xlsxwriter

def writeoutputs(stratnames, stratdict, wb, evol): 
    excelsheet = wb.add_worksheet(evol)
    k = 2
    for z in stratnames: 
        if z[1] == '0': 
            excelsheet.write(1, k, z)
            excelsheet.write(2, k, 'SPFN')
            excelsheet.write(2, k+1, 'SPFP')
            k += 2
    row = 4
    for rep in range(10):
        excelsheet.write(row, 1, "R%s" % str(rep))
        row += 1
    # write in the values 
    row, col = 3, 2
    lastval = 0
    for z in stratnames: 
        rVal = z[1]
        #print(rVal)
        if rVal != lastval: 
            lastval = rVal 
            col = 2
            row += 1
        spfn, spfp = stratdict[z]
        excelsheet.write(row, col, spfn)
        excelsheet.write(row, col+1, spfp)
        col += 2

def readoutputs(stratnames, filename): 
    stratdict = dict()
    with open(filename, "r") as f:
        lines = f.readlines()
        for i, line in enumerate(lines): 
            vals = line.split(',')
            spfn, spfp = vals[2].split('SPFN')[1][:-3], vals[3].split('SPFP')[1][:-3]
            name = stratnames[i]
            stratdict[name] = (float(spfn), float(spfp))
    return stratdict

if __name__ == '__main__': 
    # grep the SP-Score > decomp_results.txt
    filename = 'm1_10reps_compare.txt'
    exceloutname = 'm1_10reps_compare_results.xlsx'
    workbook = xlsxwriter.Workbook(exceloutname)

    evol = "M1"
    stratnames = []
    for i in range(10):  # [0]
        rep = "R%s" % str(i)
        for d in [30, 25, 20, 15, 10, 5, 1]:
            decomp = str(d)
            strats = ['stefan_fastUPP',
                    'stefan_UPP',
                        'stefan_trueUPP',
                        'stefan_UPPadjusted'
            ]

            for strat in strats: 
                name = rep + "_" + strat + "_" + evol + "_" + decomp
                stratnames.append(name)
    
    stratdict = readoutputs(stratnames, filename)
    writeoutputs(stratnames, stratdict, workbook, evol)
    workbook.close()