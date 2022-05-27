from pybedtools import BedTool
import sys
def filter_peaks(bed, log10_pval=3, log2_fc = 3):
    ''' filter peak region by pvalue and fold change for normalized peaks'''
    if log2_fc == None:
        return bed.filter(lambda x: float(x[3])>log10_pval).saveas()

    return bed.filter(lambda x: float(x[3])>log10_pval and float(x[4])>log2_fc).saveas()

if __name__=='__main__':
    inbed = sys.argv[1]
    outbed = sys.argv[2]

    filter_peaks(BedTool(inbed)).saveas(outbed)