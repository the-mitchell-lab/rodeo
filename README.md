# RODEO 2
Welcome to the new RODEO! Exhaustive documentation on how RODEO2 should be used can be found on the [RODEO website.](http://www.ripp.rodeo/index.html) The explanation below is a brief summary of how to use RODEO. Please report any bugs either via GitHub or by contacting the developers directly.

Lastly, I have one important note. IF YOU WANT TO KILL RODEO/END THE PROCESS, USE Ctrl+C! Do not use Ctrl+Z or another combination, as these are not able to be processed properly and could result in [ZOMBIE PROCESSES](https://stackoverflow.com/questions/20688982/zombie-process-vs-orphan-process). 

## Requirements
* numpy >=1.11
* SciKit Learn >=0.19
* Meme Suite >=4.12
* Biopython >=1.69
* HMMER >=3.1b2
* [RREFinder](https://github.com/Alexamk/RREFinder)
* [RADAR](https://github.com/AndreasHeger/radar)

## Installation

1. Pull the git repository down onto your computer.
2. In the hmm_dir folder of the repo, press your Pfam-A file (press TIGRFAM too if desired)
    * You can just copy the hmm files from RODEO 1.0 (Perl version) into the folder instead if you'd like.
    * If you are worried about space, you can edit the PFAM_DIR variable in the general section of the `confs/default.conf` file.

Note: If you are a MacOS user, if you are getting a "<urlopen error [SSL: CERTIFICATE_VERIFY_FAILED] certificate verify failed..." when running RODEO, please go to your applications folder and go to your Python 3 folder and run the "Install_Certificates.command" to fix this issue.

There may be things missing. Please let me know.

## Usage
I will run through a few examples and explain what they mean.
1. Basic RODEO run. 
    *  `python rodeo_main.py query`
    *  RODEO will run on query, and output will be in a folder titled `[query]_rodeo_out`
    *  `python rodeo_main.py BAL72546` for example, will run RODEO on the BAL72546 accession and output results to `BAL72546_rodeo_out`
    *  If you had a local version of AB593691 (Nucleotide contig containing BAL72546) as a .gbk file, you could run...
        * `python rodeo_main.py AB593691.gbk` 
    *  `python rodeo_main.py lassos.txt` will run on every accession in `lassos.txt` and output results to `lassos_rodeo_out`. Notice that the file extension is ignored when naming the output folder.
2. Basic RODEO with named output.
    * Say you wanted output in a particular folder, for example, `my_output`.
    * `python rodeo_main.py lassos.txt -out my_output`
    * Output will appear in `my_output_rodeo_out`.
3. RODEO with custom HMMs.
    * RODEO by default requires Pfam-A for HMM scanning. However, some RiPP hueristics make use of TIGRFAM. If you'd like to run TIGRFAM or another custom HMM, then use the `-hmm` or `--custom_hmm` flag with the path to your HMM. Note that you can input a list of HMMs as you will see below.
    * `python rodeo_main.py lassos.txt -hmm TIGRFAM.hmm MYFAVHMM.hmm` will use Pfam-A in addition to TIGRFAM and MYFAVHMM. Note that this syntax works only if the hmms are in the top level directory, as these are relative paths to rodeo_main.py.
4. Running in parallel
    * If you have a list of accessions and want to run RODEO on them in parallel, use the `-j` or `--num_cores` flag followed by the number of processes you want to spawn. Note that it doesn't make sense to spawn more processes than your computer has CPUs. It also doens't make sense to spawn 4 processes if there are only 3 queries in the input file. What's the 4th process going to do?
    * `python rodeo_main.py lassos.txt -j 4`
    * NOTE: The output may not be in the same order as the input. This is because each process will write its output as soon as possible. Let me know if this is a serious inconvenience!
5. Configuration files.
    * RODEO by default has a configuration file in `confs/` named `default.conf`. Take a look at this to get a feel for the syntax. 
    * The use of the conf file is to specify command line arguments in a file, provide colors for ORF diagrams, and make specific parameters for different types of RiPPs. If you make your own conf file in the `confs` folder, any parameter you specify will overwrite the corresponding one in the default configuration, however anything you don't specify will just use the default config value. You can supply multiple config files the same way that you supply multiple HMMs. Say you provide the flag `--conf_file confs/conf1 confs/conf2`, the parameters in conf2 take precedence over those in 1, and the parameters in conf1 take precedence over those in the default conf.
    * See the Configuration section for a more detailed description of syntax.
    * `python rodeo_main.py lassos.txt --conf_file confs/myconf.conf`

### Configuration syntax
Configurations should be placed in the conf folder. The syntax is as follows for each RiPP type. Note that the variable type (int or bool) only needs to be specified if it is not a string. If you are curious about what variable names are available, most command line arguments of RODEO are. Also note that parameters specific to genbank file mining should go in the 'general' section, as the genbank files are mined the same no matter the RiPP type, as it is a preproccessing step.
```
>[Ripp_type or 'general']
#BEGIN_VARIABLES
[variable_type] VARNAME1 VALUE1
...
...
[variable_type] VARNAMEn VALUEn
#END_VARIABLES

#BEGIN_COLORS
HMM_ANNOTATION_ID1 [color1]
...
...
HMM_ANNOTATION_IDm [colorm]
#END_COLORS
>>
```
### More Notes (some redudancy)
1. Output is currently not verbose (You will not see all debug output). For those of you who would like to see it, uncomment line 16 and comment line 17 in `rodeo_main.py`
2. Output is not in order if ran in parallel. Output should still make sense but the accessions might not be in the same order due to parallel processing.


### TODO list
1. Update sacti module search for rSAM
