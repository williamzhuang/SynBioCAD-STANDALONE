import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import re

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

    data = data.pivot_table('Values', ['Sample'], 'Gene')
    data.boxplot()

    for i,d in enumerate(data):
        y = data[d]
        x = np.random.normal(i+1, 0.04, len(y))
        plt.plot(x, y, mec='k', ms=7, marker="o", linestyle="None")

    plt.hlines(1,0,4,linestyle="--")
    plt.show()

# Methodology by Thorsten Kranz http://stackoverflow.com/questions/14399689/matplotlib-drawing-lines-between-points-ignoring-missing-data
def timePlot(data):
    data = data.sort(['Time'])
    grouped = data.groupby(['Gene', 'Time'], as_index = False).mean()
    grouped = grouped.pivot_table('Values', ['Time'], 'Gene')

    data.set_index(['Time'])

    xs = grouped.index
    for col in grouped:
        series = grouped[col]
        mask = np.isfinite(series)
        plt.plot(xs[mask], series[mask], linestyle = '-', marker = 'o')
        

    plt.show()

def askUser():
    filename = input("What is the filename of the data to be analyzed?: ")
    print("")
    print("What kind of data is it?")       # TODO: Automate determining the type of data.
    print("")
    print("\"1\" Condition")
    print("\"2\" Time")
    option = input("Your answer:")
    print("")


    if option == "1":
        data = parseExcelCondition(filename)

    elif option == "2":
        data = parseExcelTime(filename)
        timePlot(data)


askUser()