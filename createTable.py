# Writing to an excel
# sheet using Python
import xlwt, xlrd, sys, argparse, subprocess, csv
from xlwt import Workbook

print(sys.argv[1:])
parser = argparse.ArgumentParser()
parser.add_argument('--pvalList', nargs='+')
parser.add_argument('--typeList', nargs='+')
parser.add_argument('--data')
args = parser.parse_args()
print(args)

# Workbook is created
wb = Workbook()

###########################################################################
########################### Case Accuracy Table ###########################
###########################################################################

# add_sheet is used to create sheet.
sheet1 = wb.add_sheet('Case Accuracy')
sheet1.write(0,0,"Values")
# Column labels
counter = 1
for type in args.typeList:
    sheet1.write(0, counter, str(type))
    newcounter = 1

    # Grab the files corresponding to the type and place in arrays
    thefile = str(args.data)+'_'+str(type)+'_caseaccuracy.txt'
    thenum = 2*len(args.pvalList)
    COMMAND = "awk 'BEGIN{\"wc -l < "+str(thefile)+"\" | getline b}{if(NR>b-"+str(thenum)+") print}' "+str(thefile)+" > tempStats.txt"
    subprocess.call(COMMAND, shell=True)
    COMMAND = "awk 'NR%2==0 {print $NF}' tempStats.txt"
    caseaccuracy_list = map(float,subprocess.check_output(COMMAND, shell=True).splitlines())

    for SNP in args.pvalList:
        # Write the values indexed by newcounter
        sheet1.write(newcounter, counter, float(caseaccuracy_list[newcounter-1]))
        newcounter = newcounter + 1
    counter = counter + 1

# Row labels
counter = 1
for SNP in args.pvalList:
    sheet1.write(counter, 0, int(SNP))
    counter = counter + 1

###########################################################################
######################### Control Accuracy Table ##########################
###########################################################################

# add_sheet is used to create sheet.
sheet2 = wb.add_sheet('Control Accuracy')
sheet2.write(0,0,"Values")
# Column labels
counter = 1
for type in args.typeList:
    sheet2.write(0, counter, str(type))
    newcounter = 1

    # Grab the files corresponding to the type and place in arrays
    thefile = str(args.data)+'_'+str(type)+'_controlaccuracy.txt'
    thenum = 2*len(args.pvalList)
    COMMAND = "awk 'BEGIN{\"wc -l < "+str(thefile)+"\" | getline b}{if(NR>b-"+str(thenum)+") print}' "+str(thefile)+" > tempStats.txt"
    subprocess.call(COMMAND, shell=True)
    COMMAND = "awk 'NR%2==0 {print $NF}' tempStats.txt"
    controlaccuracy_list = map(float,subprocess.check_output(COMMAND, shell=True).splitlines())

    for SNP in args.pvalList:
        # Write the values indexed by newcounter
        sheet2.write(newcounter, counter, float(controlaccuracy_list[newcounter-1]))
        newcounter = newcounter + 1
    counter = counter + 1

# Row labels
counter = 1
for SNP in args.pvalList:
    sheet2.write(counter, 0, int(SNP))
    counter = counter + 1


###########################################################################
############################# F1 Score Table ##############################
###########################################################################

# add_sheet is used to create sheet.
sheet3 = wb.add_sheet('F1 Score')
sheet3.write(0,0,"Values")
# Column labels
counter = 1
for type in args.typeList:
    sheet3.write(0, counter, str(type))
    newcounter = 1

    # Grab the files corresponding to the type and place in arrays
    thefile = str(args.data)+'_'+str(type)+'_f1score.txt'
    thenum = 2*len(args.pvalList)
    COMMAND = "awk 'BEGIN{\"wc -l < "+str(thefile)+"\" | getline b}{if(NR>b-"+str(thenum)+") print}' "+str(thefile)+" > tempStats.txt"
    subprocess.call(COMMAND, shell=True)
    COMMAND = "awk 'NR%2==0 {print $NF}' tempStats.txt"
    f1score_list = map(float,subprocess.check_output(COMMAND, shell=True).splitlines())

    for SNP in args.pvalList:
        # Write the values indexed by newcounter
        sheet3.write(newcounter, counter, float(f1score_list[newcounter-1]))
        newcounter = newcounter + 1
    counter = counter + 1

# Row labels
counter = 1
for SNP in args.pvalList:
    sheet3.write(counter, 0, int(SNP))
    counter = counter + 1

###########################################################################
######################## Training Accuracy Table ##########################
###########################################################################

# add_sheet is used to create sheet.
sheet4 = wb.add_sheet('Training Accuracy')
sheet4.write(0,0,"Values")
# Column labels
counter = 1
for type in args.typeList:
    sheet4.write(0, counter, str(type))
    newcounter = 1

    # Grab the files corresponding to the type and place in arrays
    thefile = str(args.data)+'_'+str(type)+'_trainingaccuracy.txt'
    thenum = 2*len(args.pvalList)
    COMMAND = "awk 'BEGIN{\"wc -l < "+str(thefile)+"\" | getline b}{if(NR>b-"+str(thenum)+") print}' "+str(thefile)+" > tempStats.txt"
    subprocess.call(COMMAND, shell=True)
    COMMAND = "awk 'NR%2==0 {print $NF}' tempStats.txt"
    trainingaccuracy_list = map(float,subprocess.check_output(COMMAND, shell=True).splitlines())

    for SNP in args.pvalList:
        # Write the values indexed by newcounter
        sheet4.write(newcounter, counter, float(trainingaccuracy_list[newcounter-1]))
        newcounter = newcounter + 1
    counter = counter + 1

# Row labels
counter = 1
for SNP in args.pvalList:
    sheet4.write(counter, 0, int(SNP))
    counter = counter + 1

###########################################################################
######################### Testing Accuracy Table ##########################
###########################################################################

# add_sheet is used to create sheet.
sheet5 = wb.add_sheet('Testing Accuracy')
sheet5.write(0,0,"Values")
# Column labels
counter = 1
for type in args.typeList:
    sheet5.write(0, counter, str(type))
    newcounter = 1

    # Grab the files corresponding to the type and place in arrays
    thefile = str(args.data)+'_'+str(type)+'_bestaccuracy.txt'
    thenum = 2*len(args.pvalList)
    COMMAND = "awk 'BEGIN{\"wc -l < "+str(thefile)+"\" | getline b}{if(NR>b-"+str(thenum)+") print}' "+str(thefile)+" > tempStats.txt"
    subprocess.call(COMMAND, shell=True)
    COMMAND = "awk 'NR%2==0 {print $NF}' tempStats.txt"
    testingaccuracy_list = map(float,subprocess.check_output(COMMAND, shell=True).splitlines())

    for SNP in args.pvalList:
        # Write the values indexed by newcounter
        sheet5.write(newcounter, counter, float(f1score_list[newcounter-1]))
        newcounter = newcounter + 1
    counter = counter + 1

# Row labels
counter = 1
for SNP in args.pvalList:
    sheet5.write(counter, 0, int(SNP))
    counter = counter + 1

wb.save(str(args.data)+'_tableValues.xls')

# open the output csv
with open(str(args.data)+'_caseaccuracy.csv', 'wb') as myCsvfile:
    # define a writer
    wr = csv.writer(myCsvfile, delimiter="\t")

    # open the xlsx file
    myfile = xlrd.open_workbook(str(args.data)+'_tableValues.xls')
    # get a sheet
    mysheet = myfile.sheet_by_index(0)

    # write the rows
    for rownum in xrange(mysheet.nrows):
        wr.writerow(mysheet.row_values(rownum))

# open the output csv
with open(str(args.data)+'_controlaccuracy.csv', 'wb') as myCsvfile:
    # define a writer
    wr = csv.writer(myCsvfile, delimiter="\t")

    # open the xlsx file
    myfile = xlrd.open_workbook(str(args.data)+'_tableValues.xls')
    # get a sheet
    mysheet = myfile.sheet_by_index(1)

    # write the rows
    for rownum in xrange(mysheet.nrows):
        wr.writerow(mysheet.row_values(rownum))

# open the output csv
with open(str(args.data)+'_f1score.csv', 'wb') as myCsvfile:
    # define a writer
    wr = csv.writer(myCsvfile, delimiter="\t")

    # open the xlsx file
    myfile = xlrd.open_workbook(str(args.data)+'_tableValues.xls')
    # get a sheet
    mysheet = myfile.sheet_by_index(2)

    # write the rows
    for rownum in xrange(mysheet.nrows):
        wr.writerow(mysheet.row_values(rownum))

# open the output csv
with open(str(args.data)+'_trainingaccuracy.csv', 'wb') as myCsvfile:
    # define a writer
    wr = csv.writer(myCsvfile, delimiter="\t")

    # open the xlsx file
    myfile = xlrd.open_workbook(str(args.data)+'_tableValues.xls')
    # get a sheet
    mysheet = myfile.sheet_by_index(3)

    # write the rows
    for rownum in xrange(mysheet.nrows):
        wr.writerow(mysheet.row_values(rownum))

# open the output csv
with open(str(args.data)+'_testingaccuracy.csv', 'wb') as myCsvfile:
    # define a writer
    wr = csv.writer(myCsvfile, delimiter="\t")

    # open the xlsx file
    myfile = xlrd.open_workbook(str(args.data)+'_tableValues.xls')
    # get a sheet
    mysheet = myfile.sheet_by_index(4)

    # write the rows
    for rownum in xrange(mysheet.nrows):
        wr.writerow(mysheet.row_values(rownum))
