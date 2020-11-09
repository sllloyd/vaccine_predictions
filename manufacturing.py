#!/usr/bin/env python
# coding: utf-8
'''
--------------------------------------------------------------------------------
Vaccines Manufacturing Model
--------------------------------------------------------------------------------
Bryden Wood
License: MIT, see full license in LICENSE.txt
--------------------------------------------------------------------------------
Date: 2020-11-06
Authors:
    Jiabin Li
    Wynne Lim
    Stamatios Stamatiadis
    David Reader
--------------------------------------------------------------------------------
The purpose of the manufacturing model is to identify plausible scenarios 
for how long it will take to manufacture vaccines to treat the populations 
in need, and is split into two parts: Manufacturing Preparation, predicting 
when primary (drug substance) and secondary (drug product) manufacturing 
could start; and Manufacturing Capacity, predicting when enough doses are 
made to meet the needs of the target populations. The model uses output from
the Research and Development model, expert input and interviews, literature 
and manufacturing capacity survey data as inputs. The model uses Monte Carlo 
techniques to derive the dates at which the vaccine production targets will 
be met and should therefore be run many times to smooth statistical 
fluctuations.
--------------------------------------------------------------------------------
Developed for compatibility with Python 3.7.

This model functions as an importable module to R&D Vaccine Predictions model
('R&D model').

The model uses the following libraries from the default Python installation:
    random
    datetime
    json
    copy

The model uses the following additional libraries:
    numpy (tested on version 1.15.4)
    scipy (tested on version 1.1.0)
    pandas (tested on version 1.0.5)
--------------------------------------------------------------------------------
'''

# Module imports
import pandas as pd
import numpy as np
import random
import datetime
import json
import copy


# --------------------------------------------------------------------------------

def initialise(params, vaccines, manufacturing):
    '''
    This function creates and populates a set of global dataframes which can
    then be accessed from other functions. External data is loaded from JSON
    files.

    Parameters
    ----------
    None

    Returns
    ----------
    None
    '''

    # manufacturing model initialisation.

    global primary_cap
    global ratio_pri_v
    global ratio_pri_d
    global plants_categories

    global secondary_cap
    global ratio_sec
    global sec_plants_available

    global df_Categories
    global df_product

    global default
    global m_params

    global funding
    global demand
    global targetDoses

    global output_summary
    global cumulative_summary
    global df_schedule

    # load offline default data and online parameters
    default = jread(manufacturing)
    m_params = params

    primary_cap, ratio_pri_v, ratio_pri_d = getPrimaryInput()

    # initialise the product info
    df_Categories = default['Capacity Category']
    df_product = getProduct(vaccines)

    # get the number of primary plants available for each category
    plants_categories = [6, 9, 6, 0, 6]

    # initialise the schedule input
    df_schedule, funding = getScheduleInput()

    # initialise the secondary input
    secondary_cap, ratio_sec, sec_plants_available = getSecondaryInput()

    # initialise the demand data
    demand = m_params['Doses needed']
    targets = ["Group 1", "Group 2", "Group 3", "Group 4"]
    wastage = demand['Percentage of drug wastage']
    dose_perVx = demand['Number of doses per vaccine']
    # set the target numbers
    targetDoses = [demand[x] * (1 + wastage / 100) * dose_perVx for x in targets]
    targetDoses = np.cumsum(targetDoses)

    # define global variables
    output_summary = pd.DataFrame()
    cumulative_summary = pd.DataFrame()


# --------------------------------------------------------------------------------

def getSecondaryInput():
    '''
    This function extracts the relevant information pertaining to 'secondary'
    from the dictionaries: 'default' and 'm_params' and
    stores them as a separate dictionary 'secDefaultCap' and 'parameter' respectively

    Parameters
    ----------
    None

    Returns
    ----------
    secDefaultCap: a dictionary that stores the default secondary capacity of each country.

    ratio_sec: an array of ratio for low, most likely, high.

    sec_plants_available: a value for the number of secondary plants available for each category
    '''

    # secDefaultCap refers to the offline input table with the default (highest available) capacity for each country
    secDefaultCap = default['Sec Capacity']
    parameter = m_params['Sec Input']

    # get the ratio for low, most likely, high
    parameter_default = 4268.328666666668  # this is the value that is the original value we calculated

    cols = [name + ' Available Capacity (million doses/month)' for name in ['Lowest', 'Most Likely', 'Highest']]
    ratio_sec = [parameter[col] / parameter_default for col in cols]

    # get the number of secondary plants available for each category
    sec_plants_available = 0
    for key, value in secDefaultCap.items():
        if (value > 0): sec_plants_available += 1

    return secDefaultCap, ratio_sec, sec_plants_available


def getIteration(df_rnd):
    '''
    Takes the output from the R&D model and creates the iteration data required
    for the manufacturing functions. This includes joining with vaccine data in
    the database and restricting the number of successful vaccines in each
    category based on the model assumptions.

    Parameters
    ----------
    df_rnd: a dataframe of successful vaccines from the R&D model for the
        current iteration

    Returns
    ----------
    df_iteration: a dataframe of the successful vaccines with data processed
        for use in the manufacturing model
    '''

    # generate iteration table
    tryID = df_rnd['try']

    v = df_rnd['vaccines']
    df = pd.DataFrame(v, columns=['Vaccine', 1, 2, 3, 'Approval (month)'])
    df['try'] = tryID
    df = df[['try', 'Vaccine', 'Approval (month)']]

    # filter the vaccines being approved
    df_rnd = df[df['Approval (month)'] != '']
    df_rnd.reset_index(drop=True, inplace=True)

    # find platform and category
    df_iteration = df_rnd
    df_iteration = pd.merge(df_iteration, df_product, on='Vaccine', how='left')

    platforms = np.array(df_iteration['Platform'])
    categories = []

    for i in df_iteration.index.values:
        categories.append(df_Categories[platforms[i]])

    df_iteration['Category'] = categories

    # Keep the first 3 of vaccine in the same category
    df_iteration = df_iteration.groupby('Category').head(3)
    df_iteration.reset_index(drop=True, inplace=True)

    return df_iteration


# --------------------------------------------------------------------------------

def getProduct(vaccines):
    '''
    This function loads the vaccine database, selects the columns that are
    required (number, name and platform) and stores them in a dataframe.

    Parameters
    ----------
    None

    Returns
    ----------
    df_product: filtered vaccine table containing only essential columns
    '''

    # filter the vaccine data
    data = vaccines
    df_product = pd.DataFrame(columns=['Vaccine', 'Platform', 'Funding'])

    vaccines = []
    platforms = []
    funding = []

    for i in range(len(data)):
        vaccines.append(data[i]['number'])
        platforms.append(data[i]['platform'])
        funding.append(data[i]['funding_category'])

    df_product['Vaccine'] = vaccines
    df_product['Vaccine'] = df_product['Vaccine'].astype(int)
    df_product['Platform'] = platforms
    df_product['Funding'] = funding

    return df_product


# --------------------------------------------------------------------------------

def getManufacturingStartTime(df_iteration):
    '''
    Gets the manufacturing start time based on the current iteration parameters
    and the required manufacturing preparation timeline using the getSchedule
    function.

    Parameters
    ----------
    df_iteration: a dataframe of the successful vaccines with data processed
        for use in the manufacturing model

    Returns
    ----------
    None
    '''

    platforms = np.array(df_iteration['Platform'])
    approval_months = np.array(df_iteration['Approval (month)'])
    funding = np.array(df_iteration['Funding'])
    pri_starts = np.zeros(len(df_iteration))
    sec_starts = np.zeros(len(df_iteration))

    # get primary and secondary start time from scheduling (getSchedule) function
    for i in df_iteration.index.values:
        t = getSchedule(platforms[i], approval_months[i], funding[i])
        pri_starts[i] = t[0]
        sec_starts[i] = t[1]

    df_iteration['Primary Start'] = pri_starts
    df_iteration['Secondary Start'] = sec_starts
    df_iteration['Primary available'] = 0
    df_iteration['Primary assigned'] = 0

    df_iteration.drop(['Funding'], axis=1, inplace=True)


# --------------------------------------------------------------------------------
def getPrimaryInput():
    '''
    This function extracts the relevant information pertaining to 'primary'
    from the dictionaries: 'default' and 'm_params' and
    stores them as a separate dictionary 'primary_Cap' and 'parameter' respectively

    Parameters
    ----------
    None

    Returns
    ----------
    primary_cap: a dictionary that stores the default primary capacity of each country
                for each of the platform

    ratio_pri_v: a dictionary containing the ratios for the lowest, most likely and highest volume
                that is used to convert default value to online value

    ratio_pri_d: a dictionary containing the ratio for the number of doses (million doses per month)
                that is used to convert default value to online value

    '''

    # primary_cap refers to the offline input table with the default primary capacity for each country
    primary_cap = default['Pri Capacity']

    # parameter is the input table for primary capacity - from online platform
    parameter = m_params['Pri Input']

    # dfPrimary['1'], dfPrimary['2'], dfPrimary['3'], dfPrimary['4'], dfPrimary['5'] = 0, 0, 0, 0, 0

    platforms = ['DNA', 'Inactivated', 'Live-attenuated', 'Non-replicating viral vector', 'Replicating viral vector',
                 'Protein subunit', 'Other', 'RNA']
    platform_cat = ['DNA', 'Protein subunit', 'Inactivated', 'RNA']

    # dv -> default volume  dd -> default doses
    parameter_dv = [268, 1113, 779, 167728]
    parameter_dd = [112, 6488, 6488, 11596, 11596, 5798, 6355, 280]

    # ratio -> used to convert default value to online value
    ratio_pri_v = {}
    ratio_pri_d = {}
    cols = [name + ' Volume V (m3)' for name in ['Lowest', 'Most Likely', 'Highest']]

    for p in platform_cat:
        ratio_pri_v[p] = [parameter[p][col] / parameter_dv[platform_cat.index(p)] for col in cols]

    for p in platforms:
        ratio_pri_d[p] = parameter[p]['Doses for most likely volume (Million Doses per month) N_v'] / parameter_dd[
            platforms.index(p)]

    return primary_cap, ratio_pri_v, ratio_pri_d


# --------------------------------------------------------------------------------
def primary(df_iteration):
    '''
    Takes the current iteration data and allocates primary production capacity
    to each vaccine candidate based on the model assumptions.

    Parameters
    ----------
    df_iteration: a dataframe of the successful vaccines with data processed
        for use in the manufacturing model

    Returns
    ----------
    df_priAllocation: a dataframe of primary capacity allocations for the
        successful vaccines
    '''

    # sort the iteration table based on the primary start time
    df_iteration = df_iteration.sort_values(by=['Primary Start']).reset_index()

    # get capacity available for each category
    platform_cat = ['DNA', 'Protein subunit', 'Inactivated', 'RNA']
    category_list = [1, 2, 3, 5]

    platforms_bycat = {
        '1': ['DNA'],
        '2': ['Non-replicating viral vector', 'Replicating viral vector', 'Protein subunit', 'Other'],
        '3': ['Inactivated', 'Live-attenuated'],
        '5': ['RNA']
    }
    sorted_cap = {}

    primary_cap_copy = copy.deepcopy(primary_cap)

    for i in category_list:
        category_cap = primary_cap_copy[str(i)]
        example_p = platform_cat[category_list.index(i)]
        for p in platforms_bycat[str(i)]:
            ratio_cat = random.triangular(ratio_pri_v[example_p][0], ratio_pri_v[example_p][2],
                                          ratio_pri_v[example_p][1])
            for key, value in category_cap.items():
                category_cap[key][p] = category_cap[key][p] * ratio_cat

        platform_r = platforms_bycat[str(i)][0]

        sorted_temp = sorted(category_cap.items(), key=lambda x: x[1][platform_r], reverse=True)

        sorted_cap[str(i)] = sorted_temp

    # primary info initialisation

    # calculate the number of plants available per vaccine
    vx_categories = [df_iteration[df_iteration['Category'] == i].index.size for i in range(1, 6)]

    # assign no. of plants available to each vaccine product
    pri_plants_available = plants_categories[:]

    # Initialise arrays
    iter_categories = np.array(df_iteration['Category'])
    iter_pri_avail = np.array(df_iteration['Primary available'])
    iter_try = np.array(df_iteration['try'])
    iter_vaccines = np.array(df_iteration['Vaccine'])
    iter_platforms = np.array(df_iteration['Platform'])
    iter_pri_starts = np.array(df_iteration['Primary Start'])
    iter_pri_assigned = np.array(df_iteration['Primary assigned'])

    for i in range(len(df_iteration)):
        category = iter_categories[i] - 1

        if vx_categories[category] > 1:
            p = round(pri_plants_available[category] / vx_categories[category])
            iter_pri_avail[i] = p if p <= 3 else 3  # one vaccine can be assigned to maximum 3 plants
        else:
            p = pri_plants_available[category]
            iter_pri_avail[i] = p if p <= 3 else 3
        pri_plants_available[category] -= p if p <= 3 else 3
        vx_categories[category] -= 1

    df_iteration['Primary available'] = iter_pri_avail

    ### Primary product allocation
    df_priAllocation = pd.DataFrame(columns=['try', 'Country', 'Vaccine', 'Platform', 'Primary Start', 'Capacity'])

    prialloc_try = []
    prialloc_countries = []
    prialloc_vaccines = []
    prialloc_platforms = []
    prialloc_pri_starts = []
    prialloc_capacities = []

    for i in range(len(df_iteration)):
        while iter_pri_assigned[i] < iter_pri_avail[i]:
            iter_pri_assigned[i] += 1
            category = str(iter_categories[i])
            platform = iter_platforms[i]
            sorted_primary = sorted_cap[category]
            for country in sorted_primary:
                if country[1][platform] > 0:
                    prialloc_try.append(iter_try[i])
                    prialloc_countries.append(country[0])
                    prialloc_vaccines.append(iter_vaccines[i])
                    prialloc_platforms.append(iter_platforms[i])
                    prialloc_pri_starts.append(iter_pri_starts[i])
                    prialloc_capacities.append(country[1][platform] / 1000 * ratio_pri_d[platform])
                    sorted_primary.remove(country)
                    break

    df_priAllocation['try'] = prialloc_try
    df_priAllocation['Country'] = prialloc_countries
    df_priAllocation['Vaccine'] = prialloc_vaccines
    df_priAllocation['Platform'] = prialloc_platforms
    df_priAllocation['Primary Start'] = prialloc_pri_starts
    df_priAllocation['Capacity'] = prialloc_capacities

    return df_priAllocation


# --------------------------------------------------------------------------------

def secondary(df_iteration, df_priAllocation):
    '''
    Takes the current iteration data and allocates secondary production capacity
    to each vaccine candidate based on the model assumptions.

    Parameters
    ----------
    df_iteration: a dataframe of the successful vaccines with data processed
        for use in the manufacturing model
    df_priAllocation: a dataframe of primary capacity allocations for the
        successful vaccines

    Returns
    ----------
    df_secCumProduction: a dataframe of secondary capacity allocations for the
        successful vaccines, and cumulative dose count each month
    '''

    # sort the iteration table based on the secondary start time
    df_iteration = df_iteration.sort_values(by=['Secondary Start']).reset_index(drop=True)

    # get the number of plants available
    sec_plantsC = sec_plants_available
    secondary_copy = copy.deepcopy(secondary_cap)

    # get random secondary capacity for each country based on a triangular distribution based on low,most likely,high value
    ratios = random.triangular(ratio_sec[0], ratio_sec[2], ratio_sec[1])
    for key, value in secondary_copy.items():
        secondary_copy[key] = secondary_copy[key] * ratios

    # sort the secondary plants by capacity
    sort_secondary = sorted(secondary_copy.items(), key=lambda x: x[1], reverse=True)

    # get the no. of plants available for each vaccine
    v_count = df_iteration.index.size
    df_iteration['Secondary available'] = 0
    df_iteration['Secondary assigned'] = 0

    # Initialise arrays
    iter_sec_avail = np.zeros(len(df_iteration))
    iter_sec_assigned = np.zeros(len(df_iteration))
    iter_try = np.array(df_iteration['try'])
    iter_vaccines = np.array(df_iteration['Vaccine'])
    iter_platforms = np.array(df_iteration['Platform'])
    iter_sec_starts = np.array(df_iteration['Secondary Start'])

    for i in range(len(df_iteration)):
        if v_count > 1:
            p = round(sec_plantsC / v_count)
            iter_sec_avail[i] = p if p <= 3 else 3  # one vaccine can be assigned to maximum 3 plants
        else:
            p = sec_plantsC
            iter_sec_avail[i] = p if p <= 3 else 3  # one vaccine can be assigned to maximum 3 plants
        v_count -= 1
        sec_plantsC -= p

    df_iteration['Secondary available'] = iter_sec_avail

    ### secondary product allocation

    # initialisation
    df_secAllocation = pd.DataFrame(columns=['try', 'Country', 'Vaccine', 'Platform', 'Secondary Start', 'Capacity'])
    df_iteration['Secondary assigned'] = 0

    iter_sec_assigned = np.array(df_iteration['Secondary assigned'])

    secalloc_try = []
    secalloc_countries = []
    secalloc_vaccines = []
    secalloc_platforms = []
    secalloc_sec_starts = []
    secalloc_capacities = []

    # assign plants to vaccines
    for i in range(len(df_iteration)):
        while iter_sec_assigned[i] < iter_sec_avail[i]:
            iter_sec_assigned[i] += 1
            for country in sort_secondary:
                if country[1] > 0:
                    secalloc_try.append(iter_try[i])
                    secalloc_countries.append(country[0])
                    secalloc_vaccines.append(iter_vaccines[i])
                    secalloc_platforms.append(iter_platforms[i])
                    secalloc_sec_starts.append(iter_sec_starts[i])
                    secalloc_capacities.append(country[1] / 1000)
                    sort_secondary.remove(country)
                    break

    df_secAllocation['try'] = secalloc_try
    df_secAllocation['Country'] = secalloc_countries
    df_secAllocation['Vaccine'] = secalloc_vaccines
    df_secAllocation['Platform'] = secalloc_platforms
    df_secAllocation['Secondary Start'] = secalloc_sec_starts
    df_secAllocation['Capacity'] = secalloc_capacities

    # check whether primary or secondary is the bottle neck
    df_priThroughput = df_priAllocation.groupby('Vaccine')['Capacity'].sum()
    df_secThroughput = df_secAllocation.groupby('Vaccine')['Capacity'].sum()

    for i in range(len(df_secAllocation)):
        vx_ID = secalloc_vaccines[i]
        if df_secThroughput[vx_ID] > df_priThroughput[vx_ID]:
            df_secAllocation.loc[i, 'Capacity'] *= (df_priThroughput[vx_ID] / df_secThroughput[vx_ID])

    # get secondary throughput table
    df_secThroughput = df_secAllocation.copy()
    capacity_arr = np.array(df_secThroughput['Capacity'])
    start_arr = np.array(df_secThroughput['Secondary Start'])
    monthly_throughput = {}

    for j in range(1, 101):
        col_month = []
        for i in range(len(df_secThroughput)):
            if j > start_arr[i]:
                capacity_month = (j - start_arr[i]) * capacity_arr[i]
            else:
                capacity_month = 0

            col_month.append(capacity_month)

        monthly_throughput[j] = col_month

    df_monthly = pd.DataFrame(monthly_throughput)
    df_secThroughput = pd.concat([df_secThroughput, df_monthly], axis=1)

    # get cumulative production for secondary
    df_secCumProduction = df_secThroughput.copy()

    return df_secCumProduction


# --------------------------------------------------------------------------------

def getSchedule(platform, approval_month, funding_cat):
    '''
    Calculates the critical path for a vaccine's manufacturing preparation
    schedule, and returns the months that primary and secondary manufacturing
    can start.

    Parameters
    ----------
    platform: the vaccine's platform
    approval_month: the month that the vaccine reached approval from R&D

    Returns
    ----------
    m_month: the months that primary and secondary manufacturing can start
    '''

    # get the gantt table
    df_gantt = df_schedule[platform].copy()

    funding_ratio = funding[funding_cat]['Gantt duration factor*']
    # update the dependencies based on the funding category

    if funding[funding_cat]['Simultaneous tech transfer?'] == 1:  # Simultaneous Tech Transfer
        df_gantt.loc['12', 'Predecessor'] = '0'
        df_gantt.loc['22', 'Predecessor'] = '21'
    else:
        pass

    if funding[funding_cat]['Manufacturing start before approval?'] == 1:  # Manufacturing before approval
        df_gantt.loc['17', 'Predecessor'] = '17'
        df_gantt.loc['27', 'Predecessor'] = '25,27'
    else:
        pass

    # Remove rows not in use
    df_gantt['Predecessor'] = df_gantt['Predecessor'].astype(str)

    # get time for each activity
    df_gantt['Time (days)'] = df_gantt[['Type', 'Value', 'Low', 'Most Likely', 'High']].apply(getValue, axis=1)

    # apply gantt duration factor
    df_gantt['Time (days)'] = df_gantt['Time (days)'] * funding_ratio

    # Initialise columns
    df_gantt['end_date'] = datetime.date(2020, 3, 1)
    df_gantt['start_date'] = datetime.date(2020, 3, 1)

    # Calculate end time

    predecessors = np.array(df_gantt['Predecessor'])
    start_dates = np.array(df_gantt['start_date'])
    end_dates = np.array(df_gantt['end_date'])
    time_days = np.array(df_gantt['Time (days)'])

    df_gantt = df_gantt[df_gantt['Predecessor'] != '-1']

    for i in df_gantt.index.values.astype(int):
        # assign start date
        if predecessors[i] != '0':

            if str(predecessors[i]).find(',') == -1:

                pTask_index = int(df_gantt.index[df_gantt['Task ID'] == int(predecessors[i])][0])
                start_dates[i] = end_dates[pTask_index]

            else:
                temp_list = []
                for p in predecessors[i].split(','):
                    pTask_index = int(df_gantt.index[df_gantt['Task ID'] == int(p)][0])
                    temp_list.append(end_dates[pTask_index])
                start_dates[i] = max(temp_list)

            # assign end date
            if i == 5:  ## assign vaccine approval time

                date_today = datetime.date.today()
                mon = (date_today.month + approval_month) % 12

                # Assign the year and month
                if mon != 0:  # not the month of Dec
                    year = date_today.year + int((date_today.month + approval_month) / 12)
                    newmon = mon
                    nodaysinapprmonth = (datetime.date(year, newmon + 1, 1) - datetime.date(year, newmon, 1)).days

                else:  # mon=0 means that the approved month is 12
                    year = (date_today.year + int((date_today.month + approval_month) / 12)) - 1
                    newmon = 12
                    nodaysinapprmonth = (datetime.date(year + 1, 1, 1) - datetime.date(year, 12, 1)).days

                # Assign the date
                if date_today.day <= nodaysinapprmonth:
                    target_date = datetime.date(year, newmon, date_today.day)

                else:
                    target_date = datetime.date(year, newmon, nodaysinapprmonth)

                end_dates[i] = target_date
            else:
                end_dates[i] = start_dates[i] + datetime.timedelta(days=time_days[i])

        else:
            end_dates[i] = start_dates[i] + datetime.timedelta(days=time_days[i])

    m_month = []

    date_today = datetime.date.today()

    m_time = [end_dates[i] for i in [17, 27]]

    for t in m_time:
        m_month.append((t.year - date_today.year) * 12 + t.month - date_today.month)

    # apply a minimum 3 months' gap
    if m_month[0] + 3 > m_month[1]: m_month[1] = m_month[0] + 3

    return m_month


# --------------------------------------------------------------------------------


# --------------------------------------------------------------------------------

def jread(filename):
    '''
    Reads JSON files and returns the values.

    Parameters
    ----------
    filename: filename of JSON file to load

    Returns
    ----------
    values: data from the JSON file
    '''

    f = open(filename)
    values = json.load(f)
    f.close()

    return values


# --------------------------------------------------------------------------------

def getValue(values):
    '''
    Returns a single time value for an activity of the schedule based on the
    type of distribution specified.

    Parameters
    ----------
    df_schedule: dataframe of the schedule activity data
    i: row index
    x: value or distribution type

    Returns
    ----------
    computed value in days
    '''
    #### Get distribution
    return {
        'Static': values[1],
        'Triangular': random.triangular(values[2], values[4], values[3])
    }[values[0]]


# --------------------------------------------------------------------------------

def getTarget(df_iteration, df_secCumProduction):
    '''
    Finds the months at which the vaccine dose targets are hit.

    Parameters
    ----------
    df_iteration: a dataframe of the successful vaccines with data processed
        for use in the manufacturing model
    df_secCumProduction: text

    Returns
    ----------
    targetMonth: array of months when each dose target was hit
    df: vaccines for the current iteration, months the target was hit,
    '''
    targetDoses_copy = targetDoses.copy()
    targetMonth = [0, 0, 0, 0]
    # find the target month
    for i in range(1, 101):
        total = df_secCumProduction[i].sum()
        for t in range(4):
            if total > targetDoses_copy[t] and targetMonth[t] == 0:
                targetMonth[t] = i
                targetDoses_copy[t] = total
    df = df_iteration.iloc[:, :7]
    df['Target1 (month)'] = targetMonth[0]
    df['Target2 (month)'] = targetMonth[1]
    df['Target3 (month)'] = targetMonth[2]
    df['Target4 (month)'] = targetMonth[3]

    vaccine_list = df['Vaccine']
    target_dose = {}
    for k in range(4):
        t = []
        for j in range(len(df)):
            if targetMonth[k] > 0:
                t.append(df_secCumProduction.loc[df_secCumProduction['Vaccine'] == vaccine_list[j],targetMonth[k]].sum())
            else:
                t.append(0)
        target_dose['Target' + str(k+1) + ' (bn doses)'] = t

    df_target = pd.DataFrame(target_dose)
    df_target = pd.concat([df, df_target], axis=1)

    return targetMonth, df_target


# --------------------------------------------------------------------------------

def timeline(df):
    '''
    Creates the data required for the timeline bar chart.

    Parameters
    ----------
    df: dataframe of output summary data from the manufacturing model

    Returns
    ----------
    df: dataframe containing data for the timeline chart
    '''

    df = df.copy()
    df = df.iloc[:, 0:7]  # pull out columns A to F from dataframe
    df['Primary start time'] = df['Primary Start'].sub(df['Approval (month)']).astype(int)
    df['Secondary start time'] = df['Secondary Start'].sub(df['Primary Start']).astype(int)
    df['Approval (month)'] = df['Approval (month)'].astype(int)

    df = df.drop(columns=['try', 'Vaccine', 'Category', 'Primary Start', 'Secondary Start'])

    df = df.groupby(['Platform']).mean().apply(np.ceil).astype(int)

    return df


# --------------------------------------------------------------------------------

def doseBreakdown(df):
    '''
    Creates the data for the dose breakdown pie charts.

    Parameters
    ----------
    df: dataframe of output summary data from the manufacturing model

    Returns
    ----------
    target1_pie: table of data to create the pie chart for dose target 1
    target2_pie: table of data to create the pie chart for dose target 2
    target3_pie: table of data to create the pie chart for dose target 3
    target4_pie: table of data to create the pie chart for dose target 4
    '''
    #### Get dose breakdown

    df_copy = df.copy()
        
    # Create a dataframe for each Target where 0 values have been filtered out
    df_doses = [df_copy[df_copy['Target' + str(i) + ' (bn doses)'] != 0][['try','Platform','Target' + str(i) + ' (bn doses)']].reset_index(drop = True) for i in range(1,5)]
    
    # Calculate number of tries for each Target above
    tries = [i['try'].nunique() for i in df_doses]
    
    # Calculate average
    df_grouped = [(df_doses[k].groupby(['Platform'])['Target' + str(k + 1) + ' (bn doses)'].sum()) / tries[k] for k in range(len(df_doses))]

    target1, target2, target3, target4 = [i for i in df_grouped]

    # Calculate the % of each Target
    target1_percent, target2_percent, target3_percent, target4_percent = [df_grouped[i] / df_grouped[i].sum() * 100 for i in range(len(df_grouped))]
    
    df = df.groupby(['Platform'])[['Target1 (bn doses)','Target2 (bn doses)','Target3 (bn doses)','Target4 (bn doses)']].sum() / df['try'].nunique()

    target1_pie = pd.DataFrame(data = {'Target1 (bn doses)': target1, 'Target 1 %': target1_percent})
    target2_pie = pd.DataFrame(data = {'Target2 (bn doses)': target2, 'Target 2 %': target2_percent})
    target3_pie = pd.DataFrame(data = {'Target3 (bn doses)': target3, 'Target 3 %': target3_percent})
    target4_pie = pd.DataFrame(data = {'Target4 (bn doses)': target4, 'Target 4 %': target4_percent})
    
    return target1_pie, target2_pie, target3_pie, target4_pie


# --------------------------------------------------------------------------------

def getHistogram(target_no, df):
    '''
    Creates the data for the histogram chart.

    Parameters
    ----------
    target_no: vaccine dose target to generate the histogram for
    df: dataframe of output summary data from the manufacturing model

    Returns
    ----------
    return1: text
    '''
    #### Get histogram

    df = df.copy()
    # count the number of tries reach target n on a particular month
    target = pd.DataFrame(df.groupby('Target' + str(target_no) + ' (month)')['try'].nunique())

    # remove failed tries (month 0)
    if target.index[0] == 0:
        target.drop(0, inplace=True)

    # calculate the cumulative sum
    target['Cumulative %'] = np.cumsum(target.loc[:, 'try'] / sum(target['try']))

    # add referening line
    target['90%'], target['75%'], target['50%'] = (0.9, 0.75, 0.5)

    # rename the column
    target.rename(index=str, columns={"try": "Number of Runs"}, inplace=True)

    return target


# --------------------------------------------------------------------------------

def cumulativeProduction(target4_hist, cumulative_summary):
    '''
    Creates the data required for the cumulative production chart.

    Parameters
    ----------
    target4_hist: histogram data for the fourth vaccine dose target
    cumulative_summary: cumulative summary of the dose production

    Returns
    ----------
    df_cumulative: dataframe containing data for the cumulative chart
    '''
    #### Get cumulative production

    df_hist = target4_hist.copy()

    percentage = (0.1, 0.25, 0.5, 0.75, 0.9)
    month = []

    for i in percentage:
        df = df_hist[df_hist['Cumulative %'] > i]
        month.append(int(df.index[0]))

    target4 = targetDoses[3]
    df = cumulative_summary.copy()

    # Define list with column names to sum()
    lst = [i for i in range(1, 101)]

    df = df.groupby(['try'])[lst].sum()

    df_cumulative = pd.DataFrame()  # reset the dataframe
    
    for i in month:
        if i > 0:
            df1 = df[(df[i - 1] < target4) & (df[i] > target4)].head(1).T
        else:
            df1 = df[df[100]==0].head(1).T
        df_cumulative = pd.concat([df_cumulative, df1], axis=1)

    df_cumulative['Target 4'] = targetDoses[3]
    df_cumulative.reset_index(inplace=True)

    df_cumulative.columns = ['Month', '10%', '25%', '50%', '75%', '90%', 'Target 4']

    return df_cumulative


# --------------------------------------------------------------------------------

def getOutput():
    '''
    Processes the output data from the manufacturing model by calling
    functions for each individual charts required.

    Parameters
    ----------
    None

    Returns:
    ----------
    timeline_bar: table of data to create the timeline bar chart

    target1_pie: table of data to create the pie chart for dose target 1
    target2_pie: table of data to create the pie chart for dose target 2
    target3_pie: table of data to create the pie chart for dose target 3
    target4_pie: table of data to create the pie chart for dose target 4

    target1_hist: table of data to create the histogram for dose target 1
    target2_hist: table of data to create the histogram for dose target 2
    target3_hist: table of data to create the histogram for dose target 3
    target4_hist: table of data to create the histogram for dose target 4

    cum_line: table of data to create the cumulative % line chart.
    '''

        # get the timeline table
    timeline_bar = timeline(output_summary)

    # get the pie chart tables
    target1_pie, target2_pie, target3_pie, target4_pie = doseBreakdown(output_summary)

    # get the histogram tables
    target1_hist = getHistogram(1, output_summary)
    target2_hist = getHistogram(2, output_summary)
    target3_hist = getHistogram(3, output_summary)
    target4_hist = getHistogram(4, output_summary)

    # get the trendline table and highlights if targets have been met
    if target1_hist.empty:
        cum_line = {'Month': {}, '10%': {}, '25%': {}, '50%': {}, '75%': {}, '90%': {}}
        highlights = [0,0,0,0]
    else:
        cum_line = cumulativeProduction(target4_hist, cumulative_summary).to_dict()
        targets_list = [target1_hist, target2_hist, target3_hist, target4_hist]
        highlights = [int(i.index[i['Cumulative %'] >= 0.5][0]) for i in targets_list]

    return timeline_bar.to_dict(), target1_pie.to_dict(), target2_pie.to_dict(), target3_pie.to_dict(), target4_pie.to_dict() \
        , target1_hist.to_dict(), target2_hist.to_dict(), target3_hist.to_dict(), target4_hist.to_dict(), cum_line, highlights


# --------------------------------------------------------------------------------

def runTrial(trialData):
    '''
    This function takes the trial output from the R&D model and runs the main
    functions of the manufacturing model.

    Parameters
    ----------
    trialData: iteration trial data from the R&D model

    Returns
    ----------
    None
    '''

    global output_summary
    global cumulative_summary

    ## get the iteration table for this try
    df_iteration = getIteration(trialData)

    if not df_iteration.empty:
        ## append the primary and secondary manufacturing start time to the iteration table
        getManufacturingStartTime(df_iteration)

        ## get the primary allocation and throughput
        df_priAllocation = primary(df_iteration)

        ## get the secondary allocation and throughput, as well as cumulative total
        df_secCumProduction = secondary(df_iteration, df_priAllocation, )

        ## get the output for this try
        targetMonth, output_mfg = getTarget(df_iteration, df_secCumProduction)

        # merge the output table
        output_summary = output_summary.append(output_mfg, ignore_index=True)
        cumulative_summary = cumulative_summary.append(df_secCumProduction, ignore_index=True)


# --------------------------------------------------------------------------------
def getScheduleInput():
    '''
    Modifies the default Gantt chart to create a Gantt chart for each platform.

    Parameters
    ----------
    None

    Returns:
    ----------
    df_schedule:

    funding: a dictionary containing relevant information for each funding type
    '''

    # the tables below are used for value updating (static, low, most likely, high)
    gantt = pd.DataFrame(default.get('Generic Gantt')).transpose()
    parameter_v = default.get('Gantt Timelines (default)')
    parameter = m_params.get('Gantt Timelines')

    # the dictionary below is based on funding category
    funding = m_params.get('Timelines by funding criteria')

    ##############################################################################
    # Update the values of other activities to the default gantt, once only
    ranges = ['Low', 'Most Likely', 'High']

    # apply changes to all platforms
    low_d, most_likely_d, high_d = [parameter_v['All Platform'][x] for x in ranges]
    low, most_likely, high = [parameter['All Platform'][x] for x in ranges]

    low_r = low / low_d
    most_likely_r = most_likely / most_likely_d
    high_r = high / high_d

    # Convert numerical values into integers
    gantt[['Value', 'Low', 'Most Likely', 'High']] = gantt[['Value', 'Low', 'Most Likely', 'High']].astype(int)

    mask = (gantt['Activities'] != 'Scale up and and process development') & (
                gantt['Activities'] != 'Technology transfer') & (
                   gantt['Activities'] != 'DP technology transfer')

    gantt.loc[mask, 'Value'] = round(most_likely_r * gantt.loc[mask, 'Value'], 0).astype(int)
    gantt.loc[mask, 'Low'] = round(low_r * gantt.loc[mask, 'Low'], 0).astype(int)
    gantt.loc[mask, 'Most Likely'] = round(most_likely_r * gantt.loc[mask, 'Most Likely'], 0).astype(int)
    gantt.loc[mask, 'High'] = round(high_r * gantt.loc[mask, 'High'], 0).astype(int)

    gantt_updated = gantt

    del parameter['All Platform']
    ##############################################################################

    # get the iteration table for different platforms
    platforms = parameter.keys()

    # gerenate the gantt for different platform
    df_schedule = {}

    for i in platforms:
        gantt_updated_copy=gantt_updated.copy()
        low_d, most_likely_d, high_d = [parameter_v[i][x] for x in ranges]
        low, most_likely, high = [parameter[i][x] for x in ranges]

        low_r = low / low_d
        most_likely_r = most_likely / most_likely_d
        high_r = high / high_d

        mask = (gantt_updated_copy['Activities'] == 'Scale up and and process development') | (
                gantt_updated_copy['Activities'] == 'Technology transfer') | (
                           gantt_updated_copy['Activities'] == 'DP technology transfer')

        gantt_updated_copy.loc[mask, 'Low'] = round(low_r * gantt_updated_copy.loc[mask, 'Low'], 0).astype(int)
        gantt_updated_copy.loc[mask, 'Most Likely'] = round(most_likely_r * gantt_updated_copy.loc[mask, 'Most Likely'],
                                                       0).astype(int)
        gantt_updated_copy.loc[mask, 'High'] = round(high_r * gantt_updated_copy.loc[mask, 'High'], 0).astype(int)

        df_schedule[i] = gantt_updated_copy

    return df_schedule, funding

##############################################################################
