#!/usr/bin/env python
# coding: utf-8
#
# Vaccine Predictions
# Version 1.0
# models.py
#
# Simulates the COVID-19 Vaccine Pipeline using Monte Carlo methods.
# Model developed by Steve Lloyd, Anthony McDonnell, Rober Van Exan and others.
#
# Copyright 2020 Steve Lloyd
#
# Permission is hereby granted, free of charge, to any person obtaining a copy of this software
# and associated documentation files (the "Software"), to deal in the Software without restriction,
# including without limitation the rights to use, copy, modify, merge, publish, distribute,
# sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all copies or 
# substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, 
# INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
# PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE
# FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
# ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#
#===================================================================

import sys
import json
import random
import math
import datetime
import time
import traceback
import argparse
import manufacturing

#-------------------------------------------------------------------

# Set up command line arguments

parser = argparse.ArgumentParser(description='Run Vaccine Pipeline Model')
parser.add_argument('--vaccines', default='vaccines.json', help='vaccines file')
parser.add_argument('--params', default='params.json', help='parameters file')
parser.add_argument('--output', default='output.json', help='output file')
parser.add_argument('--summary', default='summary.json', help='summary file')
parser.add_argument('--csv', default='', help='trials CSV file')
args = parser.parse_args()

#-------------------------------------------------------------------

def setArgs(varg, parg, oarg, sarg, cargs=''):

	# Set arguments when using driver program
	
	args.vaccines = varg
	args.params = parg
	args.output = oarg
	args.summary = sarg
	args.csv = cargs
	
#-------------------------------------------------------------------

def runModelsMain():

	# Run the model
	
	data = {}
	Stats.init()

	try:
		
		# Read inputs from files
		
		vaccines = jread(args.vaccines)
		params = jread(args.params)

		# Check the input paramaters are sensible
		
		problems = checkParameters(params)
		if len(problems) > 0:
			string = 'Error: Incorrect Parameters:'
			for problem in problems:
				string += "\n" + problem;
			data['traceback'] = string
			raise ModelsError('Incorrect parameters')
			
		# Do some general initialisation
		
		initialise(vaccines, params)
		manufacturing.initialise(vaccines, params)

		# Set up arrays for results
								
		average = {}
		deviation = {}
		for month in range(params['months']+1):
			average[month] = {}
			deviation[month] = {}
			average[month]['Prob'] = [0] * 6
			deviation[month]['Prob'] = [0] * 6
			average[month]['Platform'] = [0] * len(params['platforms'])
			deviation[month]['Platform'] = [0] * len(params['platforms'])
			average[month]['Platform/Funding'] = []
			deviation[month]['Platform/Funding'] = []
			for platformKey in range(len(params['platforms'])):
				average[month]['Platform/Funding'].append([0] * len(params['funding_categories']))
				deviation[month]['Platform/Funding'].append([0] * len(params['funding_categories']))
			
		average['Phase Length'] = {}
		average['Phase Success'] = {}
		average['Phase Numbers'] = {}
		for platformKey in range(len(params['platforms'])):
			average['Phase Length'][platformKey] = {}
			average['Phase Success'][platformKey] = {}
			average['Phase Numbers'][platformKey] = {}
												
		trialResults = "Try ID,Vaccine,Phase I (month),Phase II (month),Phase III (month),Approval (month)\n";
		
		approvals = []
		
		# Main loop over tries
		
		tries = params['tries']
		for t in range(tries):
		
			# Initialise vaccines 
			
			for j in range(len(vaccines)):
				initialiseVaccine(vaccines[j], params)
				
			trialResult = []
			failures = {}
			successes = {}
			approved = []
			buyouts = 0
			for phase in range(len(params['phases'])):
				failures[phase] = [0] * len(params['platforms'])
				successes[phase] = [0] * len(params['platforms'])
								
			count = {}
			ended = {}
			
			# Check if vaccine has already ended a phase
			
			for j in range(len(vaccines)):
				ended[j] = {}
					
				for phase in range(len(params['phases'])):
					if vaccines[j]['end'][phase] <= 0:
						ended[j][phase] = vaccines[j]['end'][phase]
						successes[phase][vaccines[j]['platform_key']] += 1
						# This should really trigger updatePos
						# But let's assume that is taken care of in the initial POS values
				
			
			# Main loop over months
			
			firstMonth = -1
			firstMonths = [-1] * (len(params['platforms']) + 3)
			for month in range(params['months']+1):
				hasSuccessFailure = False
				hasSuccessFailures = []
				for phase in range(len(params['phases'])):
					hasSuccessFailures.append([False] * len(params['platforms']))
					
				# Limit number of Phase III trials
				
				if params['do_limit_phase3'] and month > 0: limitPhase3(inPhase3, month, vaccines, params)
												
				hasApproval = False
				
				# Main loop over vaccines (randomly)
				
				inds = range(len(vaccines))
				random.shuffle(inds)
				for i in range(len(inds)):
					j = inds[i]
					if vaccines[j]['status'] == 'Failed' or vaccines[j]['status'] == 'Approved': continue

					platformKey = vaccines[j]['platform_key']
					fundingKey = vaccines[j]['funding_key']
					Stats.addx('POS Platform ' + str(month) + ' ' + str(platformKey), vaccines[j]['overall_pos'])
					Stats.addx('POS Funding ' + str(month) + ' ' + str(fundingKey), vaccines[j]['overall_pos'])
					Stats.addx('POS Platform/Funding ' + str(month) + ' ' + str(platformKey) + ' ' + str(fundingKey), vaccines[j]['overall_pos'])

					# Skip the first month - we just count it
					if month > 0:
						# Simulate this phase for this vaccine for this month
						for phase in range(len(params['phases'])):
							# Phase ended
							if month == vaccines[j]['end'][phase]:
								Stats.addx('Phase Length ' + str(phase) + ' ' + str(platformKey), vaccines[j]['end'][phase] - vaccines[j]['start'][phase])
								Stats.addx('Phase Numbers ' + str(phase,) + ' ' + str(platformKey), 1)

								# Did it succeed?
								rand = random.uniform(0.0, 1.0)
								if rand > vaccines[j]['pos'][phase]:
									# Failed
									vaccines[j]['status'] = 'Failed'
									Stats.addx('Phase Success ' + str(phase) + ' ' + str(platformKey), 0)
									# Decide whether failure is technical or funding
									rand = random.uniform(0.0, 1.0)
									if rand < float(params['funding_tech_failure'][fundingKey]):
										# Technical rather than commercial failure
										hasSuccessFailure = True
										hasSuccessFailures[phase][platformKey] = True
										failures[phase][platformKey] += 1
									
								else:
									# Succeeded
									Stats.addx('Phase Success ' + str(phase) + ' ' + str(platformKey), 1)
									hasSuccessFailure = True
									hasSuccessFailures[phase][platformKey] = True
									successes[phase][platformKey] += 1
									ended[j][phase] = month
									if phase == 4:
										hasApproval = True
										approved.append(j)
										vaccines[j]['status'] = 'Approved'
										if params['do_update_approval']: updateApproval(approved, vaccines, params)
										if firstMonth == -1: firstMonth = month
										pl = len(params['platforms'])
										if firstMonths[platformKey] == -1: firstMonths[platformKey] = month
										if int(vaccines[j]['number']) in params['cepi'] and firstMonths[pl] == -1: firstMonths[pl] = month
										if int(vaccines[j]['number']) in params['ows'] and firstMonths[pl+1] == -1: firstMonths[pl+1] = month
										if firstMonths[pl+2] == -1: firstMonths[pl+2] = month


									if params['do_buyout_bio'] and phase == 1 and fundingKey == 4:
										# Does this Bio-tech/Academic get bought out?
										rand = random.uniform(0.0, 1.0)
										if rand < float(params['bio_buyout_fract']):
											buyouts += 1
											# Promote it to Large Pharma
											initialiseVaccine(vaccines[j], params, phase+1, month, 2)

				# End loop over vaccines
								
				# Update POS on success/failure
				if params['do_update_pos'] and hasSuccessFailure: updatePos(hasSuccessFailures, successes, failures, vaccines, params)
				
				# Count Phase III trials
				
				inPhase3 = 0
				for j in range(len(vaccines)):
					if vaccines[j]['status'] == 'Failed': continue
					if vaccines[j]['end'][3] == month: continue
					if isActive(vaccines[j], 3, month): 
						inPhase3 += 1
									
				# Initialise counters
				
				count[month] = {}
				for phase in params['phases']:
					count[month][phase] = 0
					
				count[month]['Phase I/II'] = 0
				count[month]['Phase II/III'] = 0
				count[month]['Approved'] = 0
				count[month]['Total'] = 0
				
				# Loop over vaccines counting those finished and still active

				for j in range(len(vaccines)):
					if vaccines[j]['status'] == 'Failed': continue
					if isActive(vaccines[j], 0, month): count[month]['Pre-Clinical'] += 1
					if isActive(vaccines[j], 1, month): count[month]['Phase I'] += 1
					if isActive(vaccines[j], 2, month): count[month]['Phase II'] += 1
					if isActive(vaccines[j], 3, month): count[month]['Phase III'] += 1

					if isActive(vaccines[j], 1, month) and isActive(vaccines[j], 2, month)  and (not isActive(vaccines[j], 3, month)): count[month]['Phase I/II'] += 1
					if (not isActive(vaccines[j], 1, month)) and isActive(vaccines[j], 2, month)  and isActive(vaccines[j], 3, month): count[month]['Phase II/III'] += 1
					if isActive(vaccines[j], 4, month): count[month]['Approval'] += 1

					if isFinished(vaccines[j], 4, month): count[month]['Approved'] += 1
					
					count[month]['Total'] += 1

				# End loop over vaccines
									
				# Gather statistics
				
				for phase in params['phases']:
					Stats.addx('Count ' + str(month) + ' ' + str(phase), count[month][phase])
					
				Stats.addx('Count ' + str(month) + ' Phase I/II', count[month]['Phase I/II'])
				Stats.addx('Count ' + str(month) + ' Phase II/III', count[month]['Phase II/III'])
				Stats.addx('Count ' + str(month) + ' Approved', count[month]['Approved'])
				Stats.addx('Count ' + str(month) + ' Total', count[month]['Total'])
				if count[month]['Approved'] == 0:
					Stats.addx('Count ' + str(month) + ' Prob 0', 1)
				for k in range(1, 6):
					if count[month]['Approved'] >= k:
						Stats.addx('Count ' + str(month) +  ' Prob ' + str(k), 1)
			
			# End loop over months
			
			for j in range(len(vaccines)):
				for phase in range(len(params['phases'])):
					if phase in ended[j]: 
						Stats.addx('End ' + str(j) + ' ' + str(phase), ended[j][phase])
				
			Stats.addx('Buyouts', buyouts)
			Stats.addx('First Month ' + str(firstMonth), 1)
			for i in range(len(params['platforms']) + 3):
				if firstMonths[i] != -1: Stats.addx('First Months ' + str(i), firstMonths[i])

			# Get the data for the trials csv file
			
#			if args.csv:
			trialData = {}
			for j in range(len(vaccines)):
				if len(ended[j]) == 0: continue
				entry = [t+1, int(vaccines[j]['number']), '', '', '', '']
				got = False
				for phase in range(1, len(params['phases'])):
					if phase in ended[j]:
						entry[phase+1] = ended[j][phase]
						got = True
			
				if not got: continue
				for i in range(len(entry)):
					if i > 0: trialResults += ','
					trialResults += str(entry[i])
				
				if not 'try' in trialData:
					trialData['try'] = t+1
					trialData['vaccines'] = []
					
				trialData['vaccines'].append(entry[1:])
				trialResults += "\n";
			
			approvals.append(approved)
			
			manufacturing.runTrial(trialData)
			
		# End loop over tries
				
		# Calculate averages and errors
		
		for month in range(params['months']+1):
			for phase in params['phases']:
				average[month][phase] = Stats.average('Count ' + str(month) + ' ' + str(phase))
				deviation[month][phase] = Stats.deviation('Count ' + str(month) + ' ' + str(phase))

			for platformKey in range(len(params['platforms'])):
				if Stats.count('POS Platform ' + str(month) + ' ' + str(platformKey)) > 0:
					average[month]['Platform'][platformKey] = Stats.average('POS Platform ' + str(month) + ' ' + str(platformKey))
					deviation[month]['Platform'][platformKey] = Stats.deviation('POS Platform ' + str(month) + ' ' + str(platformKey))
				for fundingKey in range(len(params['funding_categories'])):
					if Stats.count('POS Platform/Funding ' + str(month) + ' ' + str(platformKey) + ' ' + str(fundingKey)) > 0:
						average[month]['Platform/Funding'][platformKey][fundingKey] = Stats.average('POS Platform/Funding ' + str(month) + ' ' + str(platformKey) + ' ' + str(fundingKey))
				for phase in range(len(params['phases'])):
					average['Phase Length'][platformKey][phase] = Stats.average('Phase Length ' + str(phase) + ' ' + str(platformKey))
					average['Phase Success'][platformKey][phase] = Stats.average('Phase Success ' + str(phase) + ' ' + str(platformKey))
					average['Phase Numbers'][platformKey][phase] = Stats.efficiency('Phase Numbers ' + str(phase) + ' ' + str(platformKey), tries)
				
			for fundingKey in range(len(params['funding_categories'])):
				if Stats.count('POS Funding ' + str(month) + ' ' + str(fundingKey)) > 0:
					average[month]['Funding-' + str(fundingKey)] = Stats.average('POS Funding ' + str(month) + ' ' + str(fundingKey))
				
			average[month]['Phase I/II'] = Stats.average('Count ' + str(month) + ' Phase I/II')
			average[month]['Phase II/III'] = Stats.average('Count ' + str(month) + ' Phase II/III')
			average[month]['Approved'] = Stats.average('Count ' + str(month) + ' Approved')
			average[month]['Total'] = Stats.average('Count ' + str(month) + ' Total')
			
			deviation[month]['Phase I/II'] = Stats.deviation('Count ' + str(month) + ' Phase I/II')
			deviation[month]['Phase II/III'] = Stats.deviation('Count ' + str(month) + ' Phase II/III')
			deviation[month]['Approved'] = Stats.deviation('Count ' + str(month) + ' Approved')
			deviation[month]['Total'] = Stats.deviation('Count ' + str(month) +  ' Total')

			for k in range(6):
				average[month]['Prob'][k] = Stats.efficiency('Count ' + str(month) + ' Prob ' + str(k), tries)
				deviation[month]['Prob'][k] = Stats.efficiencyError('Count ' + str(month) + ' Prob ' + str(k), tries)
		
		for j in range(len(vaccines)):
			vaccines[j]['finished'] = []
			for phase in range(len(params['phases'])):
				vaccines[j]['finished'].append({})
				vaccines[j]['finished'][phase]['aveEff'] = Stats.efficiency('End ' + str(j) + ' ' + str(phase), tries)
				vaccines[j]['finished'][phase]['devEff'] = Stats.efficiencyError('End ' + str(j) + ' ' + str(phase), tries)
				vaccines[j]['finished'][phase]['aveMon'] = Stats.average('End ' + str(j) + ' ' + str(phase))
				vaccines[j]['finished'][phase]['devMon'] = Stats.deviation('End ' + str(j) + ' ' + str(phase))
		
		average['buyouts'] = Stats.average('Buyouts')
		average['first_month'] = {}
		for month in range(-1, params['months']+1):
			average['first_month'][month] = Stats.efficiency('First Month ' + str(month), tries)
		
		for j in range(len(vaccines)):
			vaccines[j]['corr'] = [0] * len(vaccines)
			
		for i in range(len(approvals)):
			for j in range(len(approvals[i])):
				for k in range(len(approvals[i])):
					vaccines[approvals[i][j]]['corr'][vaccines[approvals[i][k]]['id']] += 1
		for j in range(len(vaccines)):
			for k in range(len(vaccines)):
				vaccines[j]['corr'][k] = vaccines[j]['corr'][k]/float(params['tries'])		
		
		# Gather analysis data
		
		analysis = {}
		pl = len(params['platforms'])
		analysis['number'] = 0
		analysis['success_1'] = [0] * (pl + 3)
		analysis['success_2'] = [0] * (pl + 3)
		analysis['success_3'] = [0] * (pl + 3)
		analysis['per_run'] = [0] * (pl + 3)
		analysis['outputs'] = [0] * (pl + 3)
		analysis['months'] = [0] * (pl + 3)
		for i in range(len(approvals)):
			numb = [0] * (pl + 3)
			for k in range(len(approvals[i])):
				analysis['number'] += 1
				j = approvals[i][k]
				plat = vaccines[j]['platform_key']
				numb[plat] += 1
				numb[pl+2] += 1
				if int(vaccines[j]['number']) in params['cepi']:
					numb[pl] += 1
					
				if int(vaccines[j]['number']) in params['ows']:
					numb[pl+1] += 1

			for i in range(pl + 3):
				if numb[i] >= 1: analysis['success_1'][i] += 1
				if numb[i] >= 2: analysis['success_2'][i] += 1
				if numb[i] >= 3: analysis['success_3'][i] += 1
				analysis['per_run'][i] += numb[i]
				if numb[i] >= 1: analysis['outputs'][i] += numb[i]
				
		for i in range(pl + 3):
			if analysis['success_1'][i] > 0:
				analysis['outputs'][i] = float(analysis['outputs'][i]) / float(analysis['success_1'][i])
			analysis['success_1'][i] = float(analysis['success_1'][i]) / float(params['tries'])
			analysis['success_2'][i] = float(analysis['success_2'][i]) / float(params['tries'])
			analysis['success_3'][i] = float(analysis['success_3'][i]) / float(params['tries'])
			analysis['per_run'][i] = float(analysis['per_run'][i]) / float(params['tries'])
			analysis['months'][i] = Stats.average('First Months ' + str(i))
					
		# Gather best vaccines data
		
		bestVaccines = []
		bests = []
		cumm = 0
		for i in range(params['bestVaccinesMax']):
			num = [0] * len(vaccines)
			best = -1
			bestNum = 0
			for t in range(len(approvals)):
				tot = 0
				ok = True
				for j in bests:
					if j in approvals[t]:
						ok = False
						break
				if not ok: continue
				
				tot += 1
				for k in range(len(approvals[t])):
					j = approvals[t][k]
					num[j] += 1
					if num[j]> bestNum: 
						bestNum = num[j]
						best = j
			cumm += float(num[best])/tries
			if best > -1: bestVaccines.append({'institutes': vaccines[best]['institutes'].replace("'", ""), 'number': vaccines[best]['number'], 'id': vaccines[best]['id'], 'platform': vaccines[best]['platform'], 'funding_category': vaccines[best]['funding_category'], 'country': vaccines[best]['country'], 'percent': float(num[best])/tries, 'cumulative': cumm, 'success': vaccines[best]['finished'][4]['aveEff'], 'cepi': vaccines[best]['cepi'], 'ows': vaccines[best]['ows']})
			bests.append(best)
			
		# Sort vaccines by success
		
		def sortByRank(a, b):
			for phase in range(len(params['phases'])):
				i = len(params['phases']) - phase - 1
				comp = cmp(b['finished'][i]['aveEff'], a['finished'][i]['aveEff'])
				if comp != 0: return comp
				comp = cmp(a['finished'][i]['aveMon'], b['finished'][i]['aveMon'])
				if comp != 0: return comp
			return comp
			
		vaccines.sort(sortByRank)

		# Copy vaccines to output array
		
		vaccine_outputs = []
		for j in range(len(vaccines)):
			vaccines[j]['rank'] = j + 1
			vaccine_output = {'number': vaccines[j]['number'], 'rank': vaccines[j]['rank'],  'active': vaccines[j]['active'], 'id': vaccines[j]['id'], 'platform': vaccines[j]['platform'], 'funding_category': vaccines[j]['funding_category'], 'country': vaccines[j]['country'], 'cepi': 0, 'ows': 0}
			vaccine_output['finished'] = []
			for phase in range(len(params['phases'])):
				vaccine_output['finished'].append(vaccines[j]['finished'][phase])
			vaccine_output['institutes'] = vaccines[j]['institutes'].replace("'", "")
			vaccine_output['corr'] = vaccines[j]['corr'][:]
			vaccine_output['cepi'] = vaccines[j]['cepi']
			vaccine_output['ows'] = vaccines[j]['ows']
			vaccine_outputs.append(vaccine_output)

		# Output trials CSV
			
		if args.csv:
			f = open(args.csv, "w")
			f.write(trialResults)
			f.close()

		# Collect benchmark summary
		
		prob50 = ''
		prob90 = ''
		prob99 = ''
		for m in range(params['months']+1):
			if prob50 == '' and average[m]['Prob'][1] >= 0.5: prob50 = m
			if prob90 == '' and average[m]['Prob'][1] >= 0.9: prob90 = m
			if prob99 == '' and average[m]['Prob'][1] >= 0.99: prob99 = m
		
		prob1 = average[params['months']]['Prob'][1]
		numb = average[params['months']]['Approved']

		summ = [prob50, prob90, prob99, prob1, numb]
		jwrite(summ, args.summary)
		
		# Do some simple checks
		
		checks = []
		if params['cross_check']:
			checks = crossCheck(vaccines, params)
		
		manufacturingOutput = manufacturing.getOutput()			
		
		# Collect the output

		data['time'] = params['time_now']
		data['tries'] = params['tries']
		data['vaccines'] = vaccine_outputs
		data['average'] = average
		data['deviation'] = deviation
		data['best'] = bestVaccines
		data['analysis'] = analysis
		data['checks'] = checks
		data['manufacturing'] = manufacturingOutput
		
		status = 0

	except ModelsError:
		status = 1
		
	except:
		data['traceback'] = 'Error: ' + "\n" + traceback.format_exc()
		print 'Error:', traceback.format_exc()
		status = 1

#	print status
	
	# Write out the results

	jwrite(data, args.output)

#===================================================================

class Logger():

	# Class for collecting debug output
	
	output = ''
	
#-------------------------------------------------------------------

	@classmethod
	def getOutput(cls):
		return cls.output

#-------------------------------------------------------------------
		
	@classmethod
	def debug(cls, msg):
		cls.output += 'DEBUG: ' + msg + '<br>'
		
#===================================================================

class Stats():

	# Class for collecting statistics

#-------------------------------------------------------------------

	@classmethod
	def init(cls):
	
		cls.sumn = {}
		cls.sumx = {}
		cls.sumy = {}
		cls.sumxx = {}
		cls.sumyy = {}
		cls.sumxy = {}
	
#-------------------------------------------------------------------
	
	@classmethod
	def addx(cls, name, x):
	
		# Add a value to name
		
		if not name in cls.sumn:
			cls.sumn[name] = 0
			cls.sumx[name] = 0
			cls.sumxx[name] = 0
		
		cls.sumn[name] += 1
		cls.sumx[name] += x
		cls.sumxx[name] += x*x
		
#-------------------------------------------------------------------

	@classmethod
	def addxy(cls, name, x, y):

		if not name in cls.sumn:
			cls.sumn[name] = 0
			cls.sumx[name] = 0
			cls.sumy[name] = 0
			cls.sumxx[name] = 0
			cls.sumxy[name] = 0
			cls.sumyy[name] = 0
		
		cls.sumn[name] += 1
		cls.sumx[name] += x
		cls.sumy[name] += y
		cls.sumxx[name] += x*x
		cls.sumxy[name] += x*y
		cls.sumyy[name] += y*y
		
#-------------------------------------------------------------------

	@classmethod
	def average(cls, name):
		try:
			narg = cls.sumn[name]
			if narg == 0:
				raise StatsError('No data')

			xarg = cls.sumx[name]
						
			return float(xarg)/float(narg)

		except:
			return 0.0
		
#-------------------------------------------------------------------

	@classmethod
	def count(cls, name):
		try:
			narg = cls.sumn[name]
				
			return narg

		except:
			return 0
		
#-------------------------------------------------------------------

	@classmethod
	def countError(cls, name):
		try:
			narg = cls.sumn[name]
				
			return math.sqrt(narg)

		except:
			return 0
		
#-------------------------------------------------------------------

	@classmethod
	def deviation(cls, name):
		try:
			narg = cls.sumn[name]
			if narg < 0:
				raise StatsError('Not enough data')

			xarg = cls.sumx[name]
			xxarg = cls.sumxx[name]
							
			return math.sqrt((float(narg*xxarg) - float(xarg*xarg))/float(narg*(narg-1)))

		except:
			return 0.0
			
#-------------------------------------------------------------------

	@classmethod
	def efficiency(cls, name, total):
		try:
			if total == 0:
				raise StatsError('No data')
			
			num = cls.sumn[name]

			return float(num)/float(total)

		except:
			Logger.debug('Failed in efficiency')
			return 0
		
#-------------------------------------------------------------------

	@classmethod
	def efficiencyError(cls, name, total):

		# e = n/(n+r)
		# de^2 = (de/dn)^2*Delta_n^2 + (de/dr)^2*Delta_r^2
		# Delta_n and Delta_r are given by their upper or lower Poisson limits

		try:
			if total == 0:
				raise StatsError('No data')
			
			n = cls.sumn[name]
			r = total - n

			if n <= total/2:
				dn = cls.upper(n)
				dr = cls.lower(r)
			else:
				dn = cls.lower(n)
				dr = cls.upper(r)
			
			return math.sqrt(((r*dn)**2 + (n*dr)**2)/total**4)

		except:
			Logger.debug('Failed in efficiencyError')
			return 0
#-------------------------------------------------------------------

	@classmethod
	def upper(cls, n, S=1.0):

		# See Confidence Limits for Small Numbers of Events in Astrophysical Data
		# Neil Gehrels The Astrophysical Journal, 303:336-346, 198
		# S=1 corresponds to Gaussian statistics 1 sigma = 0.8413
		
		return S*math.sqrt(n + 3.0/4.0) + (S**2 + 3.0)/4.0
	
#-------------------------------------------------------------------

	@classmethod
	def lower(cls, n, S=1.0):

		# See Confidence Limits for Small Numbers of Events in Astrophysical Data
		# Neil Gehrels The Astrophysical Journal, 303:336-346, 198
		# S=1 corresponds to Gaussian statistics 1 sigma = 0.8413
		
		if n == 0: return 0
		return n - n*(1 - 1/(9.0*n) - S/(3.0*math.sqrt(n)))**3
		
#===================================================================

class StatsError(Exception):
	def __init__(self, value):
		self.value = value
	def __str__(self):
		return repr(self.value)
			
#===================================================================

class ModelsError(Exception):
	def __init__(self, value):
		self.value = value
	def __str__(self):
		return repr(self.value)
			
#===================================================================

def checkParameters(params):

	# Check paramaters are sensible
	
	limits = {
		'tries': ([''], 1, 50000),
		'months': ([''], 24, 240),
		'pos_factor': ([''], 1, 10),
		'timeline_factor': ([''], 1, 10),
		'best_timeline': ([0, 1, 2, 3, 4], 1, 48),
		'likely_timeline': ([0, 1, 2, 3, 4], 1, 48),
		'worst_timeline': ([0, 1, 2, 3, 4], 1, 48),
		'phase_success': ([0, 1, 2, 3, 4], 0, 1),
		'platform_pos': ([0, 1, 2, 3, 4, 5, 6, 7, 8], 0, 1),
		'platform_correlation_values': (['None', 'Low', 'Medium', 'Strong'], 0, 1),
		'timeline_factor_values': (['much_faster', 'faster', 'slightly_faster', 'normal', 'slightly_slower', 'slower', 'much_slower', 'very_much_slower'], 0, 10),
		'funding_pos': ([0, 1, 2, 3, 4], 0, 10),
		'funding_overlap_phase1_start': (['simult', 'mostly', 'phases12', 'phases23', 'consec', 'gaps'], 0, 10),
		'funding_overlap_approval_start': (['simult', 'mostly', 'phases12', 'phases23', 'consec', 'gaps'], 0, 10),
		'funding_overlap_phase_overlap': (['simult', 'mostly', 'phases12', 'phases23'], 0, 10),
		'funding_overlap_phase_gap': (['phases12', 'phases23', 'consec', 'gaps'], 0, 30),
		'funding_tech_failure': ([0, 1, 2, 3, 4], 0, 1),
		'approval_limit': ([''], 0, 100),
		'approval_pos': ([''], 0, 1),
		'approval_timeline': ([''], 1, 10),
		'phase3_slowdown_fract': ([''], 0, 1),
		'phase3_slowdown_factor': ([''], 1, 10),
		'phase3_limit': ([''], 1, 100),
		'phase3_limit_factor': ([''], 1, 10),
		'bio_buyout_fract': ([''], 0, 1)
	}
	choices = { 
		'option': ([''], ['Pessimistic', 'Normal', 'Optimistic']),
		'platform_timeline': ([0, 1, 2, 3, 4, 5, 6, 7, 8], ['much_faster', 'faster', 'slightly_faster', 'normal', 'slightly_slower', 'slower', 'much_slower', 'very_much_slower']),
		'platform_correlation': ([0, 1, 2, 3, 4, 5, 6, 7, 8], ['None', 'Low', 'Medium', 'Strong']),
		'funding_timeline': ([0, 1, 2, 3, 4], ['much_faster', 'faster', 'slightly_faster', 'normal', 'slightly_slower', 'slower', 'much_slower', 'very_much_slower']),
		'funding_overlap': ([0, 1, 2, 3, 4], ['simult', 'mostly', 'phases12', 'phases23', 'consec', 'gaps']),
		'phase3_action': ([''], ['Slow Down', 'Stop'])
	}

#	Uncommenting these should cause a Incorrect Paramaters error	
#	params['phase_success'][1] = 1.3
#	params['platform_timeline'][1] = 'blower'
#	params['funding_overlap'][1] = 'smult'
	
	problems = []
	for parameter in params:
		# Check within range
		if parameter in limits:
			name = parameter
			minVal = limits[parameter][1]
			maxVal = limits[parameter][2]
			for field in limits[parameter][0]:
				if field == '':
					value = params[parameter]
				else:
					value = params[parameter][field]
				if value < minVal:
					problems.append(parameter + ' ' + str(field) + ' value - ' + str(value) + ' < ' + str(minVal))
				if value > maxVal:
					problems.append(parameter + ' ' + str(field) + ' value - ' + str(value) + ' > ' + str(maxVal))
				
		if parameter in choices:
			name = parameter
			values = choices[parameter][1]
			for field in choices[parameter][0]:
				if field == '':
					value = params[parameter]
				else:
					value = params[parameter][field]
				if not (value in values):
					problems.append(parameter + ' ' + str(field) + ' value - ' + str(value) + ' unknown')
				
	# Special case of timelines where best <= likely <= worst
	
	for i in range(len(params['best_timeline'])):
		bestVal = params['best_timeline'][i]
		likelyVal = params['likely_timeline'][i]
		worstVal = params['worst_timeline'][i]
		if bestVal > likelyVal:
			problems.append('Best timeline for ' + params['phases'][i] + ' - ' + str(bestVal) + ' > likely timeline - ' + str(likelyVal))
		if likelyVal > worstVal:
			problems.append('Likely timeline for ' + params['phases'][i] + ' - ' + str(likelyVal) + ' > worst timeline - ' + str(worstVal))

	return problems
	
#-------------------------------------------------------------------
		
def crossCheck(vaccines, params):
	
	# Calculates the success and timelines from the parameters without any randomness
	
	checks = []
	for j in range(len(vaccines)):
		check = {}
		
		months = []
		for phase in range(len(params['phases'])):
			# Already know end
			if int(vaccines[j]['end_dates'][phase]) > 0 and int(vaccines[j]['end_dates'][phase]) < params['time_now']:
				(syear, smon) = getYearMonth(int(vaccines[j]['start_dates'][phase]))
				(eyear, emon) = getYearMonth(int(vaccines[j]['end_dates'][phase]))
				mons = emon - smon + 12*(eyear - syear)
			else:
				mons = params['best_timeline'][phase]
			
				mons += meanTimeline(params['best_timeline'][phase], params['likely_timeline'][phase], params['worst_timeline'][phase]) - params['best_timeline'][phase]

			months.append(mons)
			
		fundingOverlap = params['funding_overlap'][vaccines[j]['funding_key']]
		
		mons = 	months[0] + params['funding_overlap_phase1_start'][fundingOverlap]
		if fundingOverlap == 'simult' or 	fundingOverlap == 'mostly':
			mons += max(months[1], params['funding_overlap_phase_overlap'][fundingOverlap] + months[2], 2*params['funding_overlap_phase_overlap'][fundingOverlap] + months[3])
		elif fundingOverlap == 'phases12':
			mons += max(months[1], params['funding_overlap_phase_overlap'][fundingOverlap] + months[2]) + params['funding_overlap_phase_gap'][fundingOverlap] + months[3]
		elif fundingOverlap == 'phases23':
			mons += months[1] + params['funding_overlap_phase_gap'][fundingOverlap] + max(months[2], params['funding_overlap_phase_overlap'][fundingOverlap] + months[3])
		elif fundingOverlap == 'consec' or fundingOverlap == 'gaps':
			mons += months[1] + params['funding_overlap_phase_gap'][fundingOverlap] + months[2] + params['funding_overlap_phase_gap'][fundingOverlap] + months[3]	
		mons += params['funding_overlap_approval_start'][fundingOverlap] + months[4]

		mons = mons * params['timeline_factor_values'][params['platform_timeline'][vaccines[j]['platform_key']]]
		mons = mons * params['timeline_factor_values'][params['funding_timeline'][vaccines[j]['funding_key']]]
		check['months'] = mons
		
		check['success'] = 1.0
		for phase in range(len(params['phases'])):
			if int(vaccines[j]['end_dates'][phase]) > 0 and int(vaccines[j]['end_dates'][phase]) < params['time_now']: continue
			if phase == 3: 
				check['success'] = check['success'] * params['platform_pos'][vaccines[j]['platform_key']]
			else:
				check['success'] = check['success'] * params['phase_success'][phase]

		check['success'] = multiplyPos(check['success'], params['funding_pos'][vaccines[j]['funding_key']])

		checks.append(check)

	return checks
	
#-------------------------------------------------------------------
		
def formatDate(date):
	if date > 0:
		d = datetime.datetime.fromtimestamp(date)
		return d.strftime("%d/%m/%y")
	else:
		return ''
		
#-------------------------------------------------------------------

def getPhaseLength(vaccine, phase, params):

	# Calculates a random phase length 
	
	best_timeline = vaccine['best'][phase]
	likely_timeline = vaccine['likely'][phase]
	worst_timeline = vaccine['worst'][phase]

	# The probability function is an asymmetric triangle starting at best, 
	# peaking at likely and ending at worst.
	
	rand = 1.0
	prob = 0
	# Keep trying until the random number is below the function
	while rand > prob:
		rand = random.uniform(0.0, 1.0)
		months = random.uniform(worst_timeline, best_timeline)
		if months < likely_timeline:
			if likely_timeline - best_timeline > 0:
				prob = (months - best_timeline)/(likely_timeline - best_timeline)
			else:
				months = likely_timeline
				break
		else:
			if worst_timeline - likely_timeline > 0:
				prob = (worst_timeline - months)/(worst_timeline - likely_timeline)
			else:
				months = likely_timeline
				break

	# Randomly lengthen Phase III
		
	if phase == 3:
		rand = random.uniform(0.0, 1.0)
		if rand < params['phase3_slowdown_fract']:
			months = months * params['phase3_slowdown_factor']
	
	return int(months + 0.5)

#-------------------------------------------------------------------

def getPhaseStart(vaccine, phase, params):

	# Calculates the start of the phase according to the overlap strategy
	
	# Start Pre-clinical now
	
	if phase == 0: return 0 

	fundingOverlap = params['funding_overlap'][vaccine['funding_key']]
	
	# Start the first phase after Pre-Clinical ends
	
	if phase == 1: return vaccine['end'][phase-1] + params['funding_overlap_phase1_start'][fundingOverlap]
	
	# Start Approval after phase 3 ends
	
	if phase == 4: return vaccine['end'][phase-1] + params['funding_overlap_approval_start'][fundingOverlap]

	# Otherwise return the relevant gap and overlap values
	
	if fundingOverlap == 'simult':		# Almost simultaneous
		return vaccine['start'][phase-1] + params['funding_overlap_phase_overlap'][fundingOverlap]
	elif fundingOverlap == 'mostly':	# Mostly overlapped
		return vaccine['start'][phase-1] + params['funding_overlap_phase_overlap'][fundingOverlap]
	elif fundingOverlap == 'phases12':	# Phases 1 and 2 overlapped
		if phase == 2: return vaccine['start'][phase-1] + params['funding_overlap_phase_overlap'][fundingOverlap]
		if phase == 3: return vaccine['end'][phase-1] + params['funding_overlap_phase_gap'][fundingOverlap]
	elif fundingOverlap == 'phases23':	# Phases 2 and 3 overlapped
		if phase == 2: return vaccine['end'][phase-1] + params['funding_overlap_phase_gap'][fundingOverlap]
		if phase == 3: return vaccine['start'][phase-1] + params['funding_overlap_phase_overlap'][fundingOverlap]
	elif fundingOverlap == 'consec':	# Consecutive
		return vaccine['end'][phase-1] + params['funding_overlap_phase_gap'][fundingOverlap]
	elif fundingOverlap == 'gaps':	# Gaps between phases
		return vaccine['end'][phase-1] + params['funding_overlap_phase_gap'][fundingOverlap]
		
#-------------------------------------------------------------------

def getTimeNow():

	# Returns a timestamp for the current time
	
	return int(time.time())

#-------------------------------------------------------------------

def getYearMonth(aDate):

	# Returns the year and month of aDate
	
	(yr, mn, dy, hr, mi, sc, wd, yd, st) = time.gmtime(aDate)
	
	return (yr, mn)

#-----------------------------------------------------------------

def initialise(vaccines, params):
	
	# Initialise some parameters
	
	# Adjust the months for the platform and funding

	params['option_timeline_multiplier'] = 1.0
	params['option_pos_multiplier'] = 1.0
	if params['option'] == 'Optimistic':
		if float(params['timeline_factor']) > 0: params['option_timeline_multiplier'] = 1.0 / float(params['timeline_factor'])
		params['option_pos_multiplier'] = float(params['pos_factor'])
	elif params['option'] == 'Pessimistic':
		params['option_timeline_multiplier'] = float(params['timeline_factor'])
		if float(params['pos_factor']) > 0: params['option_pos_multiplier'] = 1.0 / float(params['pos_factor'])

	params['time_now'] = getTimeNow()
	(thisYear, thisMonth) = getYearMonth(params['time_now'])
	params['this_year'] = thisYear
	params['this_month'] = thisMonth

	# Remove dates in the future if required

	for j in range(len(vaccines)):
		vaccines[j]['original_funding_key'] = vaccines[j]['funding_key']
		if params['do_ignore_future_dates']:
			for phase in range(len(params['phases'])):
				if int(vaccines[j]['end_dates'][phase]) > params['time_now']:
					vaccines[j]['end_dates'][phase] = 0

		vaccines[j]['cepi'] = 0
		vaccines[j]['ows'] = 0
		if int(vaccines[j]['number']) in params['cepi']: vaccines[j]['cepi'] = 1
		if int(vaccines[j]['number']) in params['ows']: vaccines[j]['ows'] = 1

#-----------------------------------------------------------------

def initialiseVaccine(vaccine, params, startPhase = 0, thisMonth = 0, fundingKey = None):
	
		platformKey = vaccine['platform_key']
		if not fundingKey: fundingKey = vaccine['original_funding_key']
		vaccine['funding_key'] = fundingKey
		platformMultiplier = params['timeline_factor_values'][params['platform_timeline'][platformKey]]
		fundingMultiplier = params['timeline_factor_values'][params['funding_timeline'][fundingKey]]

		if startPhase == 0:
			vaccine['pos'] = [0] * len(params['phases'])
			vaccine['best'] = [0] * len(params['phases'])
			vaccine['likely'] = [0] * len(params['phases'])
			vaccine['worst'] = [0] * len(params['phases'])
			vaccine['start'] = [0] * len(params['phases'])
			vaccine['end'] = [0] * len(params['phases'])
			vaccine['status'] = ''
			vaccine['overall_pos'] = 1.0

		for phase in range(startPhase):
			vaccine['overall_pos'] = vaccine['overall_pos'] * vaccine['pos'][phase]
		
		for phase in range(startPhase, len(params['phases'])):
			
			# Take the PoS from the phase PoS except for Phase III which is platform dependent
			
			if phase == 3: 
				vaccine['pos'][phase] = params['platform_pos'][platformKey]
			else:
				vaccine['pos'][phase] = params['phase_success'][phase]
				
			# Recalculate phase PoS 
			# Bio-tech/Academic are assumed to be bought out after Phase I

			if params['do_buyout_bio'] and fundingKey == 4:
				vaccine['pos'][phase] = multiplyPos(vaccine['pos'][phase], params['funding_pos'][fundingKey]**(1.0/2.0))
			else:
				vaccine['pos'][phase] = multiplyPos(vaccine['pos'][phase], params['funding_pos'][fundingKey]**(1.0/5.0))
			
			vaccine['overall_pos'] = vaccine['overall_pos'] * vaccine['pos'][phase]

			# Calculate timeline parameters
			
			vaccine['best'][phase] = int(params['best_timeline'][phase] * params['option_timeline_multiplier'] * platformMultiplier * fundingMultiplier + 0.5)
			vaccine['likely'][phase] = int(params['likely_timeline'][phase] * params['option_timeline_multiplier'] * platformMultiplier * fundingMultiplier + 0.5)
			vaccine['worst'][phase] = int(params['worst_timeline'][phase] * params['option_timeline_multiplier'] * platformMultiplier * fundingMultiplier + 0.5)

			# Convert any dates into months relative to this month
	
			if int(vaccine['start_dates'][phase]) > 0:
				(year, month) = getYearMonth(int(vaccine['start_dates'][phase]))
				months = month - params['this_month'] + 12*(year - params['this_year'])
				vaccine['start'][phase] = months

			# Set know end dates (except maybe those in the future
			
			if int(vaccine['end_dates'][phase]) > 0:
				if not params['do_ignore_future_dates'] or int(vaccine['end_dates'][phase]) < params['time_now']:
					(year, month) = getYearMonth(int(vaccine['end_dates'][phase]))
					months = month - params['this_month'] + 12*(year - params['this_year'])
					vaccine['end'][phase] = months
				
			# If the end of the phase is known use that as the best date
			 
			if int(vaccine['end_dates'][phase]) > 0 and vaccine['end'][phase] > 0:
				length = vaccine['end'][phase] - vaccine['start'][phase]
				diff = length - vaccines[j]['best'][phase]
				vaccine['best'][phase] = vaccine['best'][phase] + diff
				vaccine['likely'][phase] = vaccine['likely'][phase] + diff
				vaccine['worst'][phase] = vaccine['worst'][phase] + diff

			# If the start is not set - set it according to the overlap category
								
			if int(vaccine['start_dates'][phase]) == 0 and phase > 0:
				vaccine['start'][phase] = getPhaseStart(vaccine, phase,params)
				if vaccine['start'][phase] <= thisMonth: vaccine['start'][phase] = thisMonth + 1
				
			# Calculate the end of the phase

			if int(vaccine['end_dates'][phase]) == 0:
				vaccine['end'][phase] = vaccine['start'][phase] + getPhaseLength(vaccine, phase, params)
				# Don't allow any before this month (otherwise we would have set the date already)
				if vaccine['end'][phase] <= thisMonth: vaccine['end'][phase] = thisMonth + 1
		
#-----------------------------------------------------------------

def isActive(vaccine, phase, month):

	# Says whether a vaccine is still active for a particular phase and month
	
	if vaccine['start'][phase] <= month and vaccine['end'][phase] >= month:
		return True
	else:
		return False

#-----------------------------------------------------------------

def isFinished(vaccine, phase, month):

	# Says whether a vaccine has finished a particular phase
	if vaccine['start'][phase] < month and vaccine['end'][phase] <= month:
		return True
	else:
		return False

#-------------------------------------------------------------------

def jread(filename):

	# Reads a JSON file
	
	f = open(filename)
	values = json.load(f)
	f.close()
	return values
	
#-------------------------------------------------------------------

def jwrite(values, filename, indent=4):

	# Writes a JSON file
	
	j = json.dumps(values, indent=indent)
	f = open(filename, 'w')
	print >> f, j
	f.close()

#-------------------------------------------------------------------

def limitPhase3(inPhase3, month, vaccines, params):

	# Limits or slows down the vaccines entering Phase III
	
	inds = range(len(vaccines))
	random.shuffle(inds)
	for i in range(len(inds)):
		j = inds[i]
		if vaccines[j]['status'] == 'Failed': continue
		if vaccines[j]['start'][3] != month: continue

		# If below the limit carry on
		inPhase3 += 1
		if inPhase3 <= params['phase3_limit']: continue
		
		if params['phase3_action'] == 'Stop':	
			# If we are stopping any more delay till next month
			vaccines[j]['start'][3] += 1		
			vaccines[j]['end'][3] += 1
			# Slow down approval as well		
			vaccines[j]['start'][4] += 1		
			vaccines[j]['end'][4] += 1		
		elif params['phase3_action'] == 'Slow Down':	
			# Otherwise slow it down
			months = vaccines[j]['end'][3] - vaccines[j]['start'][3]
			newMonths = int(months * params['phase3_limit_factor'] + 0.5)
			vaccines[j]['end'][3] = vaccines[j]['start'][3] + newMonths
			vaccines[j]['start'][4] += newMonths - months		
			vaccines[j]['end'][4] += newMonths - months		
			
#-------------------------------------------------------------------
		
def meanTimeline(b, l, w):

	# This is the mean, Integral(x*f(x))/Integral(f(x)), where f(x) is an asymmetric triangle 
	# starting at b (best), peaking at l (likely) and ending at w (worst) 
	# i.e. f(x) = (x-b)/(l-b) for x < l and f(x) = (w-x)/(w-l) for x > l.

	num = 0
	denom = 0
	if l > b:
		num += l**3/(3.0*(l-b)) - b**3/(3.0*(l-b)) - b*l**2/(2.0*(l-b)) + b**3/(2.0*(l-b))
		denom += (l-b)/2.0
	if w > l:
		num += w**3/(2.0*(w-l)) - w*l**2/(2.0*(w-l)) - w**3/(3.0*(w-l)) + l**3/(3.0*(w-l))
		denom += (w-l)/2.0
	
	if denom == 0:
		mean = l
	else:
		mean = num / denom
	 
	return mean

#-------------------------------------------------------------------
		
def multiplyPos(pos, multiplier):

	# Multiply the PoS by a factor
	
	# No problem if multiplier less than 1 - will converge on zero
	if multiplier <= 1.0: return pos * multiplier
	
	# Try a straightforward multiplier
	newPos = pos * multiplier
	if newPos < 1.0: return newPos

	return 1.0
	
#-------------------------------------------------------------------
		
def updateApproval(approved, vaccines, params):

	# Update approval chances of other vaccines
	
	if len(approved) < params['approval_limit']: return

	for j in range(len(vaccines)):
		if vaccines[j]['status'] == 'Failed' or vaccines[j]['status'] == 'Approved': continue
		oldPos = vaccines[j]['pos'][4]
		oldTimeline = vaccines[j]['end'][4] - vaccines[j]['start'][4]
		
		newPos = multiplyPos(oldPos, params['approval_pos'])		
		newTimeline = oldTimeline * params['approval_timeline']
	
		vaccines[j]['pos'][4] = newPos
		vaccines[j]['end'][4] = vaccines[j]['start'][4] + int(newTimeline + 0.5)
		vaccines[j]['overall_pos'] = 1.0
		for phase in range(len(params['phases'])):
			vaccines[j]['overall_pos'] = vaccines[j]['overall_pos'] * vaccines[j]['pos'][phase]
			
#-------------------------------------------------------------------
		
def updatePos(hasSuccessFailures, successes, failures, vaccines, params):

	 # Updates PoS of vaccines with the same platform on the success or failure of one

	for platformKey in range(len(params['platforms'])):
		platformCorrelation = params['platform_correlation'][platformKey]
		if platformCorrelation == 'None': continue
		correlation = params['platform_correlation_values'][platformCorrelation]
		newPos = []
		for phase in range(len(params['phases'])):
			if phase == 0: continue
			if not hasSuccessFailures[phase][platformKey]: continue
			if successes[phase][platformKey] == 0 and failures[phase][platformKey] == 0: continue
			# Only failures trigger this at Pre-Clinical, Phase I and II
			if phase < 3 and failures[phase][platformKey] == 0: continue
			if phase == 4: continue
			aggregate = successes[phase][platformKey]/float(successes[phase][platformKey] + failures[phase][platformKey])
			if phase == 3: 
				oldPos = params['platform_pos'][platformKey]
			else:
				oldPos = params['phase_success'][phase]
			
			newPos = ((aggregate - oldPos) * correlation) + oldPos
			ratio = newPos/oldPos

			# Update POS of other vaccines on success/failure of others
			for j in range(len(vaccines)):
				if vaccines[j]['status'] == 'Failed' or vaccines[j]['status'] == 'Approved': continue
				if vaccines[j]['platform_key'] != platformKey: continue
				fundingKey = vaccines[j]['funding_key']
				vaccines[j]['overall_pos'] = 1.0
				for phase2 in range(len(params['phases'])):
					if phase2 == 3: 
						oldPos = params['platform_pos'][platformKey]
					else:
						oldPos = params['phase_success'][phase2]
						
					vaccines[j]['pos'][phase2] = multiplyPos(oldPos, ratio)
					if params['do_buyout_bio'] and fundingKey == 4:
						vaccines[j]['pos'][phase2] = multiplyPos(vaccines[j]['pos'][phase2], params['funding_pos'][fundingKey]**(1.0/2.0))
					else:
						vaccines[j]['pos'][phase2] = multiplyPos(vaccines[j]['pos'][phase2], params['funding_pos'][fundingKey]**(1.0/5.0))
						
					vaccines[j]['overall_pos'] = vaccines[j]['overall_pos'] * vaccines[j]['pos'][phase2]
			
#===================================================================

if __name__ == "__main__":
	runModelsMain()

#===================================================================


	