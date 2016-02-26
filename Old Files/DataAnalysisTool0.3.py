import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import re
import seaborn as sns

def parseExcelCondition(filename):
    df = pd.read_excel(filename, index = False)
    # Placeholder columns
    df['Date'] = ""
    df['Gene'] = "a"
    df['Sample'] = "b"
    df['Condition'] = "c"
    df2 = pd.DataFrame()

    # Iterates through, extracting data from the full name. 
    for i in df.iterrows():
        _, series = i;

        fullName = series[0]
        splitName = fullName.split("_")
        date = splitName[0]
        gene = splitName[1]
        sample = splitName[2]
        condition = splitName[3]

        series.replace(series.get('Date'), date, inplace = True)
        series.replace(series.get('Gene'), gene, inplace = True)
        series.replace(series.get('Sample'), sample, inplace = True)
        series.replace(series.get('Condition'), condition, inplace = True)
        df2 = df2.append(series)

    # Reorganizes the data into something more readable, without the full name.
    df2 = df2[['Gene', 'Date', 'Sample', 'Condition', 'Values']]

    return df2

def parseExcelTime(filename):
    df = pd.read_excel(filename, index = False)
    # Placeholder columns
    df['Date'] = ""
    df['Gene'] = "a"
    df['Sample'] = "b"
    df['Time'] = 0
    df2 = pd.DataFrame()

    # Iterates through, extracting data from the full name. 
    for i in df.iterrows():
        _, series = i;

        fullName = series[0]
        splitName = fullName.split("_")
        date = splitName[0]
        gene = splitName[1]
        sample = splitName[2]
        time = splitName[3]

        # Convert time into a universal unit (min)
        timeNumber = int(re.findall(r'\d+', time)[0])
        
        if "sec" in time:
            time = timeNumber / 60
        elif "min" in time:
            time = timeNumber
        elif "hr" in time:
            time = timeNumber * 60
        elif "day" in time:
            time = timeNumber * 60 * 24
        elif "week" in time:
            time = timeNumber * 60 * 24 * 7
        elif "year" in time:
            time = timeNumber * 60 * 24 * 365

        series.replace(series.get('Date'), date, inplace = True)
        series.replace(series.get('Gene'), gene, inplace = True)
        series.replace(series.get('Sample'), sample, inplace = True)
        series.replace(series.get('Time'), time, inplace = True)


        df2 = df2.append(series)

    # Reorganizes the data into something more readable, without the full name.
    df2 = df2[['Gene', 'Date', 'Sample', 'Time', 'Values']]

    return df2


def conditionPlot(data):

    data = data.pivot_table('Values', ['Sample'], ['Gene', 'Condition'])
    answer = input("Do you want a boxplot to go with your data? (y/n): ")
    if answer == "y":
        sns.boxplot(data=data)
    sns.stripplot(data=data, size = 6, jitter = True, edgecolor = "black")
    plt.ylabel('Values')
    plt.xlabel('Gene/Condition')
    plt.show()

def timePlot(data):
    answer = input("Do you want a linear regression for your data? (y/n)")
    if answer == "n":
        timePlotScatter(data)
    else:
        timePlotLine(data)

def timePlotLine(data):
    geneNamesDict = {}
    for _, row in data.iterrows():
        geneNamesDict[row['Gene']] = 1

    data = data.pivot_table('Values', ['Sample'], ['Gene', 'Time'])
    geneList = geneNamesDict.keys()

    counter = 1
    for key in geneList:
        plt.figure(counter)     
        tempTable = data[key]
        tempTable = tempTable.T
        tempTable['Time'] = tempTable.index
        tempTable = pd.melt(tempTable, id_vars='Time')[['Time','value']]
        sns.regplot(x='Time',y='value',data=tempTable,scatter=True)
        plt.title(key)
        plt.ylabel('Values')
        plt.xlabel('Time(min)')
        counter += 1
    plt.show()

def timePlotScatter(data):
    geneNamesDict = {}
    for _, row in data.iterrows():
        geneNamesDict[row['Gene']] = 1

    data = data.pivot_table('Values', ['Sample'], ['Gene', 'Time'])
    geneList = geneNamesDict.keys()

    counter = 1
    box = input("Do you want a boxplot for each timepoint? (y/n): ")
    for key in geneList:
        plt.figure(counter)
        tempTable = data[key]
        if box == "y":
            sns.boxplot(data=tempTable)
        sns.stripplot(data=tempTable, size = 6, jitter = True, edgecolor = "black")
        plt.title(key)
        plt.ylabel('Values')
        plt.xlabel('Time(min)')
        counter += 1
    plt.show()
# Determines what kind of data is being processed.
# Returns true if the data is time data and false if condition
def isTime(filename):
    df = pd.read_excel(filename, index=False)
    for i in df.iterrows():
        _, series = i;
        splitName = series[0].split("_")
        parameter = splitName[3]
        break;

    return any(char.isdigit() for char in parameter)

def askUser():
    print("IF USING PYTHON 3, DO NOT USE QUOTES.")
    filename = input("What is the filename of the data to be analyzed?: ")

    print("")
    print("What kind of data is it?")
    print("")
    print("\"1\" Condition")
    print("\"2\" Time")
    print("Press \"Return\" to automatically determine the data type.")
    option = input("Your answer:")
    print("")
    if option == "":
        if isTime(filename):
            data = parseExcelTime(filename)
            timePlot(data)
        else:
            data = parseExcelCondition(filename)
            conditionPlot(data)

    elif option == "1":
        data = parseExcelCondition(filename)
        conditionPlot(data)
    elif option == "2":
        data = parseExcelTime(filename)
        timePlot(data)

askUser()

# Testing Files:
# conditionData.xlsx
# timeData.xlsx