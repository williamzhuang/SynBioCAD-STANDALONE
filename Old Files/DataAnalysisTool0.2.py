import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import re
from beeswarm import * 

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


def boxAndWhiskers(data):

    data = data.pivot_table('Values', ['Sample'], ['Gene', 'Condition'])
    data.boxplot()

    for i,d in enumerate(data):
        y = data[d]
        x = np.random.normal(i+1, 0.04, len(y))
        plt.plot(x, y, mec='k', ms=7, marker="o", linestyle="None")

    plt.hlines(1,0,4,linestyle="--")
    plt.ylabel('Values')
    plt.show()


def timePlot(data):
    data = data.sort(['Time'])
    # Get a list of the geneNames
    geneNamesDict = {}
    for _, row in data.iterrows():
        geneNamesDict[row['Gene']] = 1

    data = pd.pivot_table(data,values='Values',index='Time',columns=['Gene', 'Sample'])
    geneList = geneNamesDict.keys()
    # Use geneList to plot data for each individual gene.
    for key in geneList:
        tempTable = data[key]
        jitter(tempTable.index.values, tempTable.values)
        plt.title(key)
        plt.ylabel('Values')
        plt.xlabel('Time (Min)')
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

def rand_jitter(arr):
    stdev = .01*(arr.max()-arr.min())
    if len(arr.shape) == 1:
        return arr + np.random.randn(arr.shape[0]) * stdev
    else:
        return arr + np.random.randn(arr.shape[0], arr.shape[1]) * stdev

def jitter(x, y, s=20, c='b', marker='o', cmap=None, norm=None, vmin=None, vmax=None, alpha=None, linewidths=None, verts=None, hold=None, **kwargs):
    print(x.size, y.size)
    return plt.scatter(rand_jitter(x), rand_jitter(y), s=20, c='b', marker='o', cmap=None, norm=None, vmin=None, vmax=None, alpha=None, linewidths=None, verts=None, hold=None, **kwargs)

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
            boxAndWhiskers(data)

    elif option == "1":
        data = parseExcelCondition(filename)
        boxAndWhiskers(data)
    elif option == "2":
        data = parseExcelTime(filename)
        timePlot(data)

askUser()