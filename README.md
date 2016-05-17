# identify_hoogsteens

PHENIX must be sourced to use this program.

To install PHENIX :
  Go to www.phenix-online.org, look for "Download the latest official release". After downloading, follow the instructions to to install.

To source PHENIX :
  Find where you installed PHENIX and find setpaths.sh (for bash [most macs]) or setpaths.csh (for c shell). On the command line :
  $ source path/to/setpaths.sh 

To run identify_potential_purine_decoys simply provide a PDB code. The script will download the coordinates and structure factors from the PDB. YOU MUST HAVE AN INTERNET CONNECTION FOR THIS SCRIPT TGO WORK.

  $ phenix.python identify_potential_purine_decoys.py 3jxb


To keep files, including the maps, add the --keep_files flag :

  $ phenix.python identify_potential_purine_decoys.py 3jxb --keep_files

For help :
$ phenix.python identify_potential_purine_decoys.py -h


Example output :

Running this :

  $ phenix.python identify_potential_purine_decoys.py 3jxb --keep_files

```Should provide this output :
************************************ Decoys ***********************************
PDB : 3jxb    Resolution : 1.67014746437
PARAMETERS
  positive_threshold : 3
  negative_threshold : -3
  count_threshold : 4
  test_radii : 1
************ BASES ************
A,  20, ,,DA, weak_decoy_candidate
B,  40, ,,DG, super_strong_decoy_candidate
*****************************  End Decoy summary  *****************************
****************************** Decoys Candidates ******************************
PDB : 3jxb    Resolution : 1.67014746437
PARAMETERS
  positive_threshold : 3
  negative_threshold : -3
  count_threshold : 4
                    **************** BASES *****************                   

 # of super strong decoy candidates : 1 
  B,  40, ,,DG
   test 2 positive peak sample points n : 83
   test 2 positive peak sample points n : 15
   C2 negative peak sample points n : 32
   C2 negative peak sample points n : 16
   C8 negative peak sample points n : 13
   C8 negative peak sample points n : 29

 # of strong decoy candidates : 0 

 # of weak decoy candidates : 1 
  A,  20, ,,DA
   test 2 positive peak sample points n : 6
   N1 negative peak sample points n : 6
****************************** End Decoy summary ******************************
```

Comentary on the output. As you can see 3jxb has 1 super strong decoy candidates 0 strong decoy candidates and 1 weak decoy candidates. If you used --keep_files you can bring up the coordinates and thae maps. The maps are called 3jxb.pdb_2mFo-DFc_map.ccp4 and 3jxb.pdb_mFo-DFc_map.ccp4. mFo-DFc is the difference map. Upon inspection you can see that the super strong decoy candidate, DG B 40, is indeed misfit while the weak decoy candidate, DA A 20, is well fit and thus not a decoy.

