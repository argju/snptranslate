#!/usr/bin/env python

# NAME, VERSION

# from __future__ import division, print_function
import sys
import argparse
import numpy as np
import time

def readPedigree(pedfile):
    """
        Reads a pedigree from either a separate pedigree file or from the given genotype file (real=False)
        The first 3 columns have to be: sample_name, father and mother, regardless of input file
        If the number of information columns are either 1 or 2, it will assume that no parents are
        present.
    """
    ped = {}
    pedlist = []
    with open(pedfile,'r') as fin:
        count = 0
        for line in fin:
            if line.startswith('#'): continue
            l = line.strip().split()
            if len(l) == 2: name,biomat,father,mother,family,sex = l[0],l[1],'0','0','0','3'
            else: continue
            #if name not in ped:
            ped[name] = {'biomat':biomat,
                         'father':father,
                         'mother':mother,
                         'sex':sex,
                         'rank':count,
                         'children':[]}
            count += 1
            pedlist.append(name)
            #else:
            #    sys.stderr.write('%s present more than once\n' % name)
    updatePed(ped)
    return ped,pedlist

def collectPedigree(infile):
    """
        Names are not edited.
    """
    ped = {}
    pedlist = []
    with open(infile,'r') as fin:
        count = 0
        for line in fin:
            if line.startswith('#'): continue
            l = line.strip().split()
            if l[0] == 'probeset_id':
                names = l[1:]
                for count,animal in enumerate(names):
                    #animal = animal.rsplit('_',1)[0]
                    #if animal not in ped:
                    ped[animal] = {'biomat':'0',
                                 'father':'0',
                                 'mother':'0',
                                 'sex':'3',
                                 'rank':count,
                                 'children':[]}
                    pedlist.append(animal)
                    #else:
                    #    sys.stderr.write('%s present more than once\n' % animal)
    return ped,pedlist

def updatePed(ped):
    # Assign children to parents and set sex of parents
    for n in ped:
        father,mother = ped[n]['father'],ped[n]['mother']
        if father in ped:
            #ped[father]['sex'] = '1'
            ped[father]['children'].append(n)
        if mother in ped:
            #ped[mother]['sex'] = '0'
            ped[mother]['children'].append(n)

def readMarkers(markerfile):
    """
        Read a marker file, format is Plink map-file with 4 additional columns:
        col 5+6: The two SNP alleles as 1/2/3/4
        col 7: The ID of the SNP on the array
        col 8: If the SNP should be included in the output (defualt: yes = 1) [optional]
    """
    mark = {}
    marklist = []
    with open(markerfile,'r') as fin:
        count = 0
        for line in fin:
            if line.startswith('#'): continue
            l = line.strip().split()
            keep = '1'
            if len(l) == 0: continue
            if len(l) > 6:
                chrom,name,distance,position,a1,a2,probeid = l[0:7]
                if len(l) == 8:
                    keep = l[7]
                else:
                    keep = 1
            else:
                raise Exception('Map file needs column 5, 6, 7 and 8 to be the two alleles, probe_id and if the SNP should be in the output')
            if name not in mark:
                try:
                    pos = int(position)
                except ValueError:
                    pos = 0
                mark[name] = {'chrom':chrom,
                              'pos':pos,
                              'alleles': a1+a2,
                              'probe_id': probeid,
                              'rank':count,
                              'keep':keep}
                count += 1
                mark[probeid] = name
                marklist.append(name)
        #print("Read %d markers" % count)
    return mark,marklist

def importFile(infile,ped,pedlist,mark,marklist):
    """
        Reads in all genotypes for the markers provided in the markerlist
        Converts from 0/1/2 to 1/2/3/4
        Names are not edited.
    """
    gen = np.zeros((len(pedlist),len(marklist)))
    with open(infile,'r') as fin:
        count = 0
        for line in fin:
            if line.startswith('#'): continue
            l = line.strip().split()
            if l[0] == 'probeset_id':
                names = l[1:]
                continue
            prb,genos = l[0],l[1:]
            try:
                marker = mark[prb]
            except KeyError:
                continue
            m = mark[marker]['alleles']
            icol = mark[marker]['rank']
            for ind,name in enumerate(names):
                #name = name.rsplit('_',1)[0]
                try:
                    irow = ped[name]['rank']
                except KeyError:
                    continue
                el = genos[ind]
                if el == '0':
                    gen[irow,icol] = int(m[0]+m[0])
                elif el == '1':
                    gen[irow,icol] = int(m[0]+m[1])
                elif el == '2':
                    gen[irow,icol] = int(m[1]+m[1])
    return gen

#***************************************************************************************
#*
#*
#***************************************************************************************

def writeSequenom(outfile,genos,ped,mark,pedlist,marklist,callrate=False):
    """ Writes a Sequenom matrix-report file """
    def trans(a):
        l = []
        for s in a:
            if s in ['A','1']: l.append('A')
            elif s in ['C','2']: l.append('C')
            elif s in ['G','3']: l.append('G')
            elif s in ['T','4']: l.append('T')
            else: l.append('0')
        return l

    def fix(s):
        if s == 'D': return 'DEL'
        if s == 'I': return 'INS'
        return s

    def checkDiagnose(data):
        out = '\t'
        # alpha-casein
        a1,a2 = data['CSN1S1_192']
        if a1 == 'A': out += 'A'
        elif a1 == 'G' : out += 'B'
        else: out += '0'
        if a2 == 'A': out += 'A'
        elif a2 == 'G' : out += 'B'
        else: out += '0'
        out += '\t'
        # beta-casein
        a1,a2 = data['CSN2_67']
        if a1 == 'T': out += 'A1'
        elif a1 == 'G': out += 'A2'
        else: out += '0'
        if a2 == 'T': out += '/A1'
        elif a2 == 'G': out += '/A2'
        else: out += '0'
        out += '\t'
        # kappa-casein
        a1,a2 = data['CSN3_136']
        if a1 == 'C': out += 'A'
        elif a1 == 'T': out += 'B'
        else: out += '0'
        if a2 == 'C': out += 'A'
        elif a2 == 'T': out += 'B'
        else: out += '0'
        out += '\t'
        # Shrimp-flavour
        a1,a2 = data['FMO3']
        if a1+a2 in ['CC']: out += 'HN'
        elif a1+a2 in ['CT','TC']: out += 'HB'
        elif a1+a2 in ['TT']: out += 'HM'
        else: out += '0'
        out += '\t'
        # Polled
        a1,a2 = data['BovineHD0100000698']
        b1,b2 = data['BovineHD0100000798']
        c1,c2 = data['BovineHD0100000821']
        d1,d2 = data['ss86313348']
        cat1,cat2 = 0,0
        if a1 == a2 == 'A': cat1 += 1
        elif a1 != a2: cat2 += 1
        if b1 == b2 == 'T': cat1 += 1
        elif b1 != b2: cat2 += 1
        if c1 == c2 == 'T': cat1 += 1
        elif c1 != c2: cat2 += 1
        if d1 == d2 == 'T': cat1 += 1
        elif d1 != d2: cat2 += 1
        s = a1+a2+b1+b2+c1+c2+d1+d2
        if [a1,a2] == ['C','C'] or [b1,b2] == ['C','C'] or [c1,c2] == ['G','G'] or [d1,d2] == ['C','C']: out += 'HH'
        elif cat1 == 4: out += 'KK'
        else:
            if s.count('A') > 4: out += '0'
            else: out += '0'
        out += '\t'
        # Deletion BTA12 (not in any way elegant)
        a1,a2 = data['ARS-BFGL-NGS-21160']
        b1,b2 = data['ARS-BFGL-NGS-104500']
        c1,c2 = data['Hapmap30016-BTA-149504']
        d1,d2 = data['BTA-17536']
        e1,e2 = data['Hapmap43970-BTA-92461']
        f1,f2 = data['ARS-BFGL-NGS-119219']
        g1,g2 = data['ARS-BFGL-NGS-80012']
        h1,h2 = data['ARS-BFGL-NGS-107895']
        if a1 != a2 or b1 != b2 or c1 != c2 or d1 != d2 or e1 != e2 or f1 != f2 or g1 != g2 or h1 != h2:
            out += '2'
        else: out += '1'
        return out

    import datetime
    dato = datetime.date.today().strftime("%d.%m.%Y")
    sep = '\t'
    fout = open(outfile,'w')
    totgenotype,maf,m = {},{},{}
    allele3,once = False,False
    totsamples = len(pedlist)
    for marker in marklist:
        maf[marker] = {}
        m[marker] = []
    extra = 'Prodnr_Indnr\tBiomat_ID\tDato'
    if 'CSN1S1_192' in m: extra += '\talpha-kasein'
    if 'CSN2_67' in m: extra += '\tbeta-kasein'
    if 'CSN3_136' in m: extra += '\tkappa-kasein'
    if 'FMO3' in m: extra += '\tfish-smell'
    if 'BovineHD0100000698' in m: extra += '\tPolled'
    if 'ARS-BFGL-NGS-21160' in m: extra += '\tBTA12_deletion'
    if len(marklist) > 0:
        #for ma in marklist:
        #    fout.write('\t%s' % ma)
        fout.write('%s\n' % extra)
    for animal in pedlist:
        irow = ped[animal]['rank']
        diagnose = {}
        cr = [0,0]
        fout.write('%s\t%s\t%s' % (animal,ped[animal]['biomat'],dato))
        for marker in marklist:
            icol = mark[marker]['rank']
            a = genos[irow,icol]
            if a == 0:
                a = '00'
            else:
                a = '%.0f' % a
            a1,a2 = trans(a)
            #fout.write(sep+fix(a1)+fix(a2)) # Printing
            cr[1] += 1
            if a1 != '0' and a2 != '0': cr[0] += 1
            if marker in ['BovineHD0100000698','BovineHD0100000798','BovineHD0100000821',
                          'CSN1S1_192','CSN2_106','CSN2_67','CSN3_136','CSN3_148',
                          'FMO3','LGB','ss86313348','ARS-BFGL-NGS-21160','ARS-BFGL-NGS-104500',
                          'Hapmap30016-BTA-149504','BTA-17536','Hapmap43970-BTA-92461',
                          'ARS-BFGL-NGS-119219','ARS-BFGL-NGS-80012','ARS-BFGL-NGS-107895']:
                diagnose[marker]= [a1,a2]
            if a1 != '0': maf[marker][(a1+a2)] = maf[marker].get((a1+a2),0) + 1
            if a1 != '0':
                totgenotype[marker] = totgenotype.get(marker,0) + 1
                if a1 not in m[marker]: m[marker].append(a1)
                if a2 not in m[marker]: m[marker].append(a2)
            if len(m[marker]) == 3:
                allele3 = True
                print m[marker]
        if len(diagnose) > 0:
            fout.write(checkDiagnose(diagnose))
        fout.write('\n')
    out1,out2,out3,out4 = 'Total_genotyped','NoCall','Total_no._Samples','Total_efficiency'
    if allele3:
        out5,out6,out7,out8,out9,out10 = '01:01','01:02','02:02','01:01','01:02','02:02'
        out14,out15,out16,out17,out18,out19= '01:03','02:03','03:03','01:03','02:03','03:03'
        out11,out12,out13 = 'Frequency_01:01','Frequency_01:02','Frequency_02:02'
        out20,out21,out22 = 'Frequency_01:03','Frequency_02:03','Frequency_03:03'
        for marker in marklist:
            if len(m[marker]) > 0: m1 = m[marker][0]
            else: m1 = '0'
            if len(m[marker]) > 1: m2 = m[marker][1]
            else: m2 = '0'
            if len(m[marker]) > 2: m3 = m[marker][2]
            else: m3 = '0'
            try:
                out1 += sep+str(totgenotype[marker])
                out2 += sep+str(totsamples-totgenotype[marker])
                out3 += sep+str(totsamples)
            except KeyError:
                out1 += sep+'0'
                out2 += sep+str(totsamples)
                out3 += sep+str(totsamples)
            try: out4 += sep+str(totgenotype[marker]/float(totsamples))
            except (ZeroDivisionError,KeyError): out4 += sep+'0'
            # Writing alleles
            out5 += sep+fix(m1)+fix(m1)
            out6 += sep+fix(m1)+fix(m2)
            out7 += sep+fix(m2)+fix(m2)
            out14 += sep+fix(m1)+fix(m3)
            out15 += sep+fix(m2)+fix(m3)
            out16 += sep+fix(m3)+fix(m3)
            maa = maf[marker].get((m1+m1),0)
            mab = maf[marker].get((m1+m2),0)+maf[marker].get((m2+m1),0)
            mbb = maf[marker].get((m2+m2),0)
            mac = maf[marker].get((m1+m3),0)+maf[marker].get((m3+m1),0)
            mbc = maf[marker].get((m2+m3),0)+maf[marker].get((m3+m2),0)
            mcc = maf[marker].get((m3+m3),0)
            out8 += sep+str(maa)
            out9 += sep+str(mab)
            out10 += sep+str(mbb)
            out17 += sep+str(mac)
            out18 += sep+str(mbc)
            out19 += sep+str(mcc)
            if maa+mab+mbb+mac+mbc+mcc > 0:
                out11 += sep+str(float(maa)/(maa+mab+mbb+mac+mbc+mcc))
                out12 += sep+str(float(mab)/(maa+mab+mbb+mac+mbc+mcc))
                out13 += sep+str(float(mbb)/(maa+mab+mbb+mac+mbc+mcc))
                out20 += sep+str(float(mac)/(maa+mab+mbb+mac+mbc+mcc))
                out21 += sep+str(float(mbc)/(maa+mab+mbb+mac+mbc+mcc))
                out22 += sep+str(float(mcc)/(maa+mab+mbb+mac+mbc+mcc))
            else:
                out11 += sep+'0'
                out12 += sep+'0'
                out13 += sep+'0'
                out20 += sep+'0'
                out21 += sep+'0'
                out22 += sep+'0'
        #fout.write('\n\n'+out5+'\n'+out6+'\n'+out7+'\n'+
        #           out14+'\n'+out15+'\n'+out16+'\n\n'+
        #           out8+'\n'+out9+'\n'+out10+'\n'+
        #           out17+'\n'+out18+'\n'+out19+'\n\n')
        #fout.write(out1+'\n'+out2+'\n\n'+out3+'\n\n'+out4+'\n\n\n\n'+
        #           out11+'\n'+out12+'\n'+out13+'\n'+
        #           out20+'\n'+out21+'\n'+out22+'\n')
    else:
        out5,out6,out7,out8,out9,out10 = '01:01','01:02','02:02','01:01','01:02','02:02'
        out11,out12,out13 = 'Frequency_01:01','Frequency_01:02','Frequency_02:02'
        for marker in marklist:
            if len(m[marker]) > 0: m1 = m[marker][0]
            else: m1 = '0'
            if len(m[marker]) > 1: m2 = m[marker][1]
            else: m2 = '0'
            try:
                out1 += sep+str(totgenotype[marker])
                out2 += sep+str(totsamples-totgenotype[marker])
                out3 += sep+str(totsamples)
            except KeyError:
                out1 += sep+'0'
                out2 += sep+str(totsamples)
                out3 += sep+str(totsamples)
            try: out4 += sep+str(totgenotype[marker]/float(totsamples))
            except (ZeroDivisionError,KeyError): out4 += sep+'0'
            # Printing alleles
            out5 += sep+fix(m1)+fix(m1)
            out6 += sep+fix(m1)+fix(m2)
            out7 += sep+fix(m2)+fix(m2)
            maa = maf[marker].get((m1+m1),0)
            mab = maf[marker].get((m1+m2),0)+maf[marker].get((m2+m1),0)
            mbb = maf[marker].get((m2+m2),0)
            out8 += sep+str(maa)
            out9 += sep+str(mab)
            out10 += sep+str(mbb)
            if maa+mab+mbb > 0:
                out11 += sep+str(float(maa)/(maa+mab+mbb))
                out12 += sep+str(float(mab)/(maa+mab+mbb))
                out13 += sep+str(float(mbb)/(maa+mab+mbb))
            else:
                out11 += sep+'0'
                out12 += sep+'0'
                out13 += sep+'0'
        #fout.write('\n\n'+out5+'\n'+out6+'\n'+out7+'\n\n'+
        #            out8+'\n'+out9+'\n'+out10+'\n\n'+
        #            out1+'\n'+out2+'\n\n'+out3+'\n\n'+out4+'\n\n\n\n'+
        #            out11+'\n'+out12+'\n'+out13+'\n')
    fout.close()

def writeGS(outfile,genos,ped,mark,pedlist,marklist):
    """
    Writes Geno
    """
    with open(outfile,'w') as fout:
        fout.write('#\t%s\n' % '\t'.join([m for m in marklist if mark[m]['keep'] != '0']))
        for animal in pedlist:
            ra = ped[animal]['rank']
            father = ped[animal]['father']
            mother = ped[animal]['mother']
            fout.write('%s\t%s\t%s' % (animal,father,mother))
            for i,m in enumerate(marklist):
                try:
                    if mark[m]['keep'] == '0': continue
                    rm = mark[m]['rank']
                    a = genos[ra,rm]
                except KeyError:
                    continue
                if a == 0:
                    a = '00'
                else:
                    a = '%.0f' % a
                try:
                    fout.write('\t%s\t%s' % (a[0],a[1]))
                except:
                    sys.stderr.write('ERROR: Unexpected allele: [%s,%s,%s]\n' % (animal,m,a))
                    sys.exit(1)
            fout.write('\n')

def main():
    parser = argparse.ArgumentParser(description='Read and process a Bovine Affy 50k calls-file.')
    parser.add_argument('infile',help='Axiom calls-file')
    parser.add_argument('-o','--outfile',help='Genotype file')
    parser.add_argument('-r','--reportfile',help='Report file')
    parser.add_argument('-m','--markers',help='List of markers')
    parser.add_argument('-p','--pedigree',help='List of samples')
    parser.add_argument('-v','--verbose',action='store_true',help='Prints runtime info')
    args = parser.parse_args()
    try:
        if args.markers:
            mark,marklist = readMarkers(args.markers)
        else:
            sys.stderr.write('ERROR: Missing marker information\n')
            sys.exit(1)
    except IOError:
        sys.stderr.write('ERROR: Could not open marker file\n')
        sys.exit(1)
    try:
        if args.pedigree:
            ped,pedlist = readPedigree(args.pedigree)
        else:
            ped,pedlist = collectPedigree(args.infile)
    except IOError:
        sys.stderr.write('Could not find sample information\n')
        sys.exit(1)
    genos = importFile(args.infile,ped,pedlist,mark,marklist)
    seqmarklist = ['BovineHD0100000698','BovineHD0100000798','BovineHD0100000821',
                          'CSN1S1_192','CSN2_106','CSN2_67','CSN3_136','CSN3_148',
                          'FMO3','LGB','ss86313348','ARS-BFGL-NGS-21160','ARS-BFGL-NGS-104500',
                          'Hapmap30016-BTA-149504','BTA-17536','Hapmap43970-BTA-92461',
                          'ARS-BFGL-NGS-119219','ARS-BFGL-NGS-80012','ARS-BFGL-NGS-107895']
    writeSequenom(args.reportfile,genos,ped,mark,pedlist,seqmarklist)
    writeGS(args.outfile,genos,ped,mark,pedlist,marklist)

if __name__ == '__main__':
    #t = time.time()
    main()
    #sys.stdout.write('Time spent: %.3f\n' % (time.time()-t))
