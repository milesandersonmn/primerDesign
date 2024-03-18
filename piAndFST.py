import allel
import numpy as np
import pandas
import csv
samples = ["LA2750_N02", "LA2749_NO2", "LA2932_N02", "LA4109_NO1", "LA4108_T30", "LA4336_NO1", "LA4107_N01", "LA4106_NO2", "LA4338_NO2", "LA4339_NO1", "LA2884_NO5", "LA4329_NO2", "LA4330_N02", "LA4332_X1", "LA4117A_N05", "LA4118_N02", "LA4119_V16", "LA4334_NO1", "LA4335_NO1", "LA2880_N17"]

with open('all_genes.bed') as bedFile:
	for gene in bedFile:
		try:
			region = gene.split()
			Chr = region[0]
			Start = int(region[1]) + 1
			End = int(region[2]) + 1
			Coordinate = Chr + ":" + str(Start) + "-" + str(End)
			callset = allel.read_vcf('Schil_capture_variants_AllChroms.decomposedVariants.filtered.maxMissing.sorted.reheader.vcf.gz', fields='*', region=Coordinate)
			pos = allel.SortedIndex(callset['variants/POS'])
			gt_array = allel.GenotypeArray(callset['calldata/GT'])
			alleleCountArray = gt_array.count_alleles()
			pi = allel.sequence_diversity(pos, alleleCountArray)
			callset = allel.read_vcf('Schil_capture_variants_AllChroms.decomposedVariants.filtered.maxMissing.sorted.reheader.vcf.gz', fields='*', region=Coordinate, samples=samples)
			pos = allel.SortedIndex(callset['variants/POS'])
			gt_array = allel.GenotypeArray(callset['calldata/GT'])
			subpops = [[0,4,7,8,9,10,13,14,15,19],[1,2,3,5,6,11,12,16,17,18]]
			acCoast = gt_array.count_alleles(subpop=subpops[0])
			acHighland = gt_array.count_alleles(subpop=subpops[1])
			num, den = allel.hudson_fst(acCoast, acHighland)
			fst = np.sum(num) / np.sum(den)
			Coordinate
			pi
			fst
			with open('piAndFST.csv', 'a', newline='') as csvfile:
				writer = csv.writer(csvfile, delimiter='\t')
				row = [Chr,Start,End,pi,fst]
				print(row)
				writer.writerow(row)
						
		except:
			pass
