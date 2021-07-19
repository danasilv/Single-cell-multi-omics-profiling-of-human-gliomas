----------------Preprocessing---------------
CBS_DNACopy.R: Runs DNACopy's implementation of a circular binary segmentation algorithm to detect copy number change
--RNA
	prepFilesForInferCNV.py: Takes scRNA counts and preprocesses them for inferCNV
	inferCNVBySample.R:  Runs inferCNV
	getInferCNVOutputBySample.R: takes the inferCNV output and splits out copy number

--Methylation
	CNVDataGenerator.py: Takes the raw CpGs counts by window and outputs copy number estimates
	convert_CNV_files_for_R_DNACopy.py: Takes the output from CNVDataGenerator.py and gets them in the right format for DNACopy

