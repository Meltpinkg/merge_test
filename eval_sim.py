import sys
import argparse
import logging
import time

import vcf

def phase_bnd(alt, chr, pos):
	# print(alt)
	# print(str(alt[0]))
	if alt[0] == ']':
		form = ']]N'
		chr2 = alt.split(':')[0][1:]
		pos2 = int(alt.split(':')[1][:-2])
	elif alt[0] == '[':
		form = '[[N'
		chr2 = alt.split(':')[0][1:]
		pos2 = int(alt.split(':')[1][:-2])
	else:
		# print(type(alt))
		if alt[1] == ']':
			form = 'N]]'
			chr2 = alt.split(':')[0][2:]
			pos2 = int(alt.split(':')[1][:-1])
		else:
			form = 'N[['
			chr2 = alt.split(':')[0][2:]
			pos2 = int(alt.split(':')[1][:-1])

	try:
		if int(chr) <= int(chr2):
			if form == 'N[[':
				form = ']]N'
			if form == ']]N':
				form = 'N[['
			return form, chr, pos, chr2, pos2
		else:
			return form, chr2, pos2, chr, pos
	except:
		return form, chr, pos, chr2, pos2

def load_callset(path):
	callset = dict()
	variant_call_file = open(path, 'r')
	vcf_reader = vcf.Reader(variant_call_file)
	for record in vcf_reader:
		chr = record.CHROM
		start = record.POS
		filter = record.FILTER
		try:
			svtype = record.INFO['SVTYPE'][:3]
		except:
			svtype = str(record.ALT[0])[1:4]
		try:
			svlen = record.INFO['SVLEN']
			try:
				if len(svlen):
					svlen = svlen[0]
			except:
				pass
		except:
			svlen = 0
		try:
			end = record.INFO['END']
		except:
			end = start
		try:
			genotype = record.samples[0]['GT']
		except:
			continue


		if svtype not in callset:
			callset[svtype] = dict()

		if chr not in callset[svtype]:
			callset[svtype][chr] = list()

		if svtype != "BND":
			if svlen == 0:
				svlen = end - start + 1

			callset[svtype][chr].append([start, end, abs(svlen), genotype, filter, 0])
		else:
			form, chr, start, chr2, pos2 = phase_bnd(str(record.ALT[0]), chr, start)
			if chr not in callset[svtype]:
				callset[svtype][chr] = list()
			callset[svtype][chr].append([start, chr2, pos2, genotype, filter, form, 0])

	variant_call_file.close()
	return callset

def eval(call, ans, bias, offect, genotype):
	for svtype in call:
		if svtype not in ans:
			continue
		else:
			for chr in call[svtype]:
				if chr not in ans[svtype]:
					continue
				else:
					for call_record in call[svtype][chr]:
						for ans_record in ans[svtype][chr]:
							# print(call_record, ans_record)
							if svtype == 'INS':
								if abs(call_record[0] - ans_record[0]) <= offect and float(min(call_record[2], ans_record[2]) / max(call_record[2], ans_record[2])) >= bias:
									call_record[-1] = 1
									ans_record[-1] = 1
									if call_record[3] in genotype[chr]:
										call_record[-1] = 2
										ans_record[-1] = 2
									break
							else:
								if max(call_record[0] - offect, ans_record[0]) <= min(call_record[1] + offect, ans_record[1]) and float(min(call_record[2], ans_record[2]) / max(call_record[2], ans_record[2])) >= bias:
									call_record[-1] = 1
									ans_record[-1] = 1
									if call_record[3] in genotype[chr]:
										call_record[-1] = 2
										ans_record[-1] = 2
									break

	# for svtype in call:
	# 	if svtype not in ans:
	# 		# continue
	# 		if svtype == 'INS':
	# 			for i in call[svtype]:
	# 				for key in ans:
	# 					for j in ans[key]:
	# 						if i[0] == j[0]:
	# 							if abs(i[1] - j[1]) <= offect and float(min(i[3],j[3])/max(i[3],j[3])) >= bias:
	# 								i[-1] = 1
	# 								j[3+opt] = 1
	# 								if i[4] == genotype[j[0]]:
	# 									i[-1] = 2
	# 									j[3+opt] = 2
	# 	else:
	# 		for i in call[svtype]:
	# 			for j in ans[svtype]:
	# 				if i[0] != j[0]:
	# 					continue
	# 				else:
	# 					if svtype in ['INS']:
	# 						if abs(i[1] - j[1]) <= offect and float(min(i[3],j[2])/max(i[3],j[2])) >= bias:
	# 							j[2+opt] = 1
	# 							i[-1] = 1
	# 							if i[4] == genotype[j[0]]:
	# 								j[2+opt] = 2
	# 								i[-1] = 2
	# 					elif svtype == 'BND':
	# 						if i[2] != j[2]:
	# 							continue
	# 						else:
	# 							if i[4] == j[4]:
	# 							# if 1:
	# 								if abs(i[1]-j[1]) <= offect and abs(i[3]-j[3]) <= offect:
	# 									i[-1] = 1
	# 									j[4+opt] = 1
	# 									if i[5] == genotype[j[0]] or i[5] == genotype[j[2]]:
	# 										i[-1] = 2
	# 										j[4+opt] = 2
	# 					else:
	# 						if max(i[1]-offect, j[1]) <= min(i[2]+offect, j[2]) and float(min(i[3],j[3])/max(i[3],j[3])) >= bias:
	# 							j[3+opt] = 1
	# 							i[-1] = 1
	# 							if i[4] == genotype[j[0]]:
	# 								j[3+opt] = 2
	# 								i[-1] = 2
						# genotype

def init_stalen():
	stalen = {'0-99': [0, 0, 0, 0],
			'100-499': [0, 0, 0, 0],
			'500-999': [0, 0, 0, 0],
			'1k-5k': [0, 0, 0, 0],
			'5k-10k': [0, 0, 0, 0],
			'other': [0, 0, 0, 0]}
	# base	TN 	call 	TP
	return stalen

def statistics(call, ans, res):
	for svtype in ans:
	# for svtype in ["DEL"]:
		tp = 0
		total_record = 0
		th_len = init_stalen()

		if svtype not in call:
			continue

		for chr in call[svtype]:
			for ele in call[svtype][chr]:
				total_record += 1
				if abs(ele[2]) < 100:
					th_len['0-99'][2] += 1
				elif abs(abs(ele[2])) < 500:
					th_len['100-499'][2] += 1
				elif abs(abs(ele[2])) < 1000:
					th_len['500-999'][2] += 1
				elif abs(abs(ele[2])) < 5000:
					th_len['1k-5k'][2] += 1
				elif abs(abs(ele[2])) < 10000:
					th_len['5k-10k'][2] += 1
				else:
					th_len['other'][2] += 1

				if ele[-1] >= res:
					tp += 1
					if abs(abs(ele[2])) < 100:
						th_len['0-99'][3] += 1
					elif abs(abs(ele[2])) < 500:
						th_len['100-499'][3] += 1
					elif abs(abs(ele[2])) < 1000:
						th_len['500-999'][3] += 1
					elif abs(abs(ele[2])) < 5000:
						th_len['1k-5k'][3] += 1
					elif abs(abs(ele[2])) < 10000:
						th_len['5k-10k'][3] += 1
					else:
						th_len['other'][3] += 1

		logging.info('TP-%d of %s:\t%d\t%d'%(res, svtype, tp, total_record))

		fn = 0
		total_base = 0

		for chr in ans[svtype]:
			for ele in ans[svtype][chr]:
				total_base += 1
				if abs(abs(ele[2])) < 100:
					th_len['0-99'][0] += 1
				elif abs(abs(ele[2])) < 500:
					th_len['100-499'][0] += 1
				elif abs(abs(ele[2])) < 1000:
					th_len['500-999'][0] += 1
				elif abs(abs(ele[2])) < 5000:
					th_len['1k-5k'][0] += 1
				elif abs(abs(ele[2])) < 10000:
					th_len['5k-10k'][0] += 1
				else:
					th_len['other'][0] += 1

				if ele[-1] >= res:
					fn += 1
					if abs(abs(ele[2])) < 100:
						th_len['0-99'][1] += 1
					elif abs(abs(ele[2])) < 500:
						th_len['100-499'][1] += 1
					elif abs(abs(ele[2])) < 1000:
						th_len['500-999'][1] += 1
					elif abs(abs(ele[2])) < 5000:
						th_len['1k-5k'][1] += 1
					elif abs(abs(ele[2])) < 10000:
						th_len['5k-10k'][1] += 1
					else:
						th_len['other'][1] += 1
					
		logging.info('TN-%d of %s:\t%d\t%d'%(res, svtype, fn, total_base))
		# print("%d\t%d\t%d\t%d"%(tp, total_record, fn, total_base))
		if res == 1:
			logging.info("++++++++  %s  ++++++++"%svtype)
		else:
			logging.info("++++++++  %s-GT  ++++++++"%svtype)
		logging.info('\t\tAcc\tRecall\tF1')
		try:
			logging.info('%s:\t%.4f\t%.4f\t%.4f\t(%d\t%d\t%d\t%d)'%('ALL',
						float(tp/total_record),
						float(fn/total_base),
						float(2*float(tp/total_record)*float(fn/total_base)/(float(tp/total_record)+float(fn/total_base))),
						tp,
						total_record,
						fn,
						total_base))
		except:
			logging.info('%s:\t%.4f\t%.4f\t%.4f\t(%d\t%d\t%d\t%d)'%('ALL',
						float(0),
						float(0),
						float(0),
						tp,
						total_record,
						fn,
						total_base))
		for key in th_len:
			# print(th_len[key])
			try:
				acc = float(th_len[key][3]/th_len[key][2])
			except:
				acc = 0
			try:
				recall = float(th_len[key][1]/th_len[key][0])
			except:
				recall = 0
			try:
				f1 = float(2*float(th_len[key][3]/th_len[key][2])*float(th_len[key][1]/th_len[key][0])/(float(th_len[key][3]/th_len[key][2])+th_len[key][1]/th_len[key][0]))
			except:
				f1 = 0
			logging.info('%s:\t%.4f\t%.4f\t%.4f\t(%d\t%d\t%d\t%d)'%(key,
																	acc,
																	recall,
																	f1,
																	th_len[key][3],
																	th_len[key][2],
																	th_len[key][1],
																	th_len[key][0]))
			# print("%d\t%d\t%d\t%d"%(th_len[key][3], th_len[key][2], th_len[key][1], th_len[key][0]))

typetrans = {'insertion':'INS', 
			'deletion':'DEL', 
			'inversion':'INV',
			'tandem duplication':'DUP',
			'reciprocal translocation':'BND'
			}

def load_ans(path):
	ansbed = dict()
	file = open(path, 'r')
	for line in file:
		seq = line.strip('\n').split('\t')
		chr = seq[0]
		svtype = typetrans[seq[3]]
		start = int(seq[1])
		end = int(seq[2])

		if svtype not in ansbed:
			ansbed[svtype] = dict()
			ansbed[svtype][chr] = list()
		else:
			if chr not in ansbed[svtype]:
				ansbed[svtype][chr] = list()

		if svtype == 'DEL':
			ansbed[svtype][chr].append([start, end, end - start + 1, 0])

		if svtype == 'DUP':
			ansbed[svtype][chr].append([start, end, end - start + 1, 0])

		if svtype == 'INV':
			ansbed[svtype][chr].append([start, end, end - start + 1, 0])

		if svtype == 'INS':
			ansbed[svtype][chr].append([start, end, len(seq[4]), 0])

		'''
		if svtype in ['INS']:
			ansbed[svtype].append([chr, start, len(seq[4]), 0, 0, 0, 0])
		elif svtype in ['BND']:
			chr2 = seq[4].split(':')[1]
			start2 = int(seq[4].split(':')[2])
			strand1 = seq[4].split(':')[3]
			strand2 = seq[4].split(':')[4]

			if strand1[0] == 'f':
				if strand2[0] == 'f':
					ansbed[svtype].append([chr, start, chr2, start2, "N[[", 0, 0, 0, 0])
					ansbed[svtype].append([chr, start, chr2, start2, "]]N", 0, 0, 0, 0])
					ansbed[svtype].append([chr, end, chr2, start2+end-start, "]]N", 0, 0, 0, 0])
					ansbed[svtype].append([chr, end, chr2, start2+end-start, "N[[", 0, 0, 0, 0])
				else:
					ansbed[svtype].append([chr, start, chr2, start2, "N[[", 0, 0, 0, 0])
					ansbed[svtype].append([chr, start, chr2, start2+end-start, "[[N", 0, 0, 0, 0])
					ansbed[svtype].append([chr, end, chr2, start2, "N]]", 0, 0, 0, 0])
					ansbed[svtype].append([chr, end, chr2, start2+end-start, "]]N", 0, 0, 0, 0])
			else:
				if strand2[0] == 'f':
					ansbed[svtype].append([chr, start, chr2, start2+end-start, "N]]", 0, 0, 0, 0])
					ansbed[svtype].append([chr, start, chr2, start2, "]]N", 0, 0, 0, 0])
					ansbed[svtype].append([chr, end, chr2, start2, "[[N", 0, 0, 0, 0])
					ansbed[svtype].append([chr, end, chr2, start2+end-start, "N[[", 0, 0, 0, 0])
				else:
					ansbed[svtype].append([chr, start, chr2, start2+end-start, "N]]", 0, 0, 0, 0])
					ansbed[svtype].append([chr, end, chr2, start2, "N]]", 0, 0, 0, 0])
					ansbed[svtype].append([chr, end, chr2, start2, "[[N", 0, 0, 0, 0])
					ansbed[svtype].append([chr, start, chr2, start2+end-start, "[[N", 0, 0, 0, 0])
		else:
			ansbed[svtype].append([chr, start, end, end-start+1, 0, 0, 0, 0])
		'''

	file.close()
	return ansbed

def load_gt(path):
	GT = dict()
	file = open(path, 'r')
	for line in file:
		seq = line.strip('\n').split('\t')
		chr = seq[0]
		if chr not in GT:
			GT[chr] = ''
		if float(seq[-1]) > 80.0:
			GT[chr] = ['1/1']
		elif 80.0 >= float(seq[-1]) > 20.0:
			GT[chr] = ['0/1']
		else:
			GT[chr] = ['0/0', './.']
	return GT

def main_ctrl(args):

	# load ground truth set
	ans = load_ans(args.ans)
	genotype = load_gt(args.gt)
		
	call = load_callset(args.call)
	# for key in call:
	# 	for chr in call[key]:
	# 		print(key, chr, len(call[key][chr]))
	# 		for i in call[key][chr]:
	# 			print(i)
	# logging.info("The number of calls within abnormal SV type in callset:")
	# for key in abcall:
	# 	logging.info("<call-%s>\t%d."%(key, abcall[key]))
	logging.info("Evaluation on call callsets...")
	eval(call, ans, args.bias, args.offect, genotype)
	statistics(call, ans, 1)
	statistics(call, ans, 2)


def main(argv):
	args = parseArgs(argv)
	setupLogging(False)
	# print args
	starttime = time.time()
	main_ctrl(args)
	logging.info("Finished in %0.2f seconds."%(time.time() - starttime))

USAGE="""\
	Evaluate SV callset generated by simulations.
	Author: Tao Jiang
	Email: tjiang@hit.edu.cn
"""

def parseArgs(argv):
	parser = argparse.ArgumentParser(prog="eval_sim.py", description=USAGE, formatter_class=argparse.RawDescriptionHelpFormatter)
	parser.add_argument("ans", type=str, help="Ground truth of simulations.")
	parser.add_argument("gt", type=str, help="Genotype in each chromosome.")
	parser.add_argument('call', type=str, help = "SV callsets")
	parser.add_argument('-b', '--bias', help = "Bias of overlaping.[%(default)s]", default = 0.7, type = float)
	parser.add_argument('-o', '--offect', help = "Offect of translocation overlaping.[%(default)s]", default = 1000, type = int)
	args = parser.parse_args(argv)
	return args

def setupLogging(debug=False):
	logLevel = logging.DEBUG if debug else logging.INFO
	logFormat = "%(asctime)s [%(levelname)s] %(message)s"
	logging.basicConfig( stream=sys.stderr, level=logLevel, format=logFormat )
	logging.info("Running %s" % " ".join(sys.argv))

if __name__ == '__main__':
	main(sys.argv[1:])
