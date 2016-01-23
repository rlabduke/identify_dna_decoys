import os,sys
import argparse
import re
import diff_density_utils

def get_sample_points_kin(bases,positive_threshold,negative_threshold) :
  kin = ''
  for base in bases :
    #print base.density_sample_points.sample_points;sys.exit()
    num = base.res_num.strip()
    cid = base.chain_id
    #if (num == '8' and cid == 'C'):# or (num == '8' and cid == 'C') :
    #  print base.get_string_id() + ':' + str(base.xyz.test_1)
    #  print base.get_string_id() + ':' + str(base.xyz.test_2) 
    #  kin += base.get_density_sample_points_kin_group(positive_threshold,
                                                      # negative_threshold)
    kin += base.get_density_sample_points_kin_group(positive_threshold,
                                                      negative_threshold)

  return kin

def run(args) :
  # defaults
  pos_thresh = 3
  neg_thresh = -3
  n_points = 4
  # process args
  desc = "This program identifies potential purine decoys in DNA. In this "
  desc+= "context, a purine decoy is one that may need to switch from anti to "
  desc+= "cis or vise vesa. A potential decoy is identified by having "
  desc+= "difference density peaks in expected locations around the purine. "
  desc+= "If you're thinking that a desciption of those locations should be "
  desc+= "included right here, you're correct - go shame Bradley about not "
  desc+= "completing this program description properly." 
  parser = argparse.ArgumentParser(description=desc)
  parser.add_argument('pdb_id', help='A PDB code')
  hs = 'The %s sigma value where difference density points at or %s '
  hs+= 'than that value will be flagged as %s peaks. Default = %i'
  parser.add_argument('-p','--positive_threshold',
                      help=hs % ('positive','greater','positive',pos_thresh),
                      type=int,default=pos_thresh)
  parser.add_argument('-n','--negative_threshold',
                      help=hs % ('negative','less','negative',neg_thresh),
                      type=int,default=neg_thresh)
  hs = 'The smallest number of adjacent difference density sample points '
  hs+= 'to be considered a difference density peak'
  parser.add_argument('-s','--point_n',help=hs,type=int,default=n_points)
  hs = 'A file path for output, defaults to stdout'
  parser.add_argument('-o','--out_file_name',help=hs)
  hs = 'Write out sample points in kinemage format'
  parser.add_argument('-k','--write_kin',help=hs)
  parser.add_argument('--keep_files',help='Keep generated pdb files')
  args = parser.parse_args()
  # does the provided pdb_id match the expexted PDB code format?
  es = 'The PDB code provided doesn\'t match the expexted PDB code format.'
  assert re.match('^[0-9]{1}[a-z0-9]{3}$',args.pdb_id.lower()), es
  # ensure that the tresholds are the proper signs
  assert args.positive_threshold > 0,'positive_threshold must be > 0'
  assert args.negative_threshold < 0,'negative_threshold must be < 0'
  assert args.point_n > 2, 'ponit_n must be > 2'
  # end process args

  # set log
  if args.out_file_name :
    log = open(args.out_file_name,'w')
  else : log = sys.stdout

  ddab = diff_density_utils.DiffDensityAroundBases(args.pdb_id,
                                positive_threshold = args.positive_threshold,
                                negative_threshold = args.negative_threshold,
                                count_threshold = args.point_n)

  if args.write_kin :
    fn = 'sample_points_%s.kin' % pdb_id
    fle = open(fn,'w')
    fle.write(get_sample_points_kin(ddab.dna_bases, positive_threshold,
                                              negative_threshold))
    fle.close()

  # write summary of diffence poits to stdout
  ddab.set_potential_decoys()
  if not args.keep_files : ddab.clean_up_files()
  ddab.write_potential_decoy_summary()

  if args.out_file_name : log.close()

  print 'success'

if __name__ == '__main__' :
  run(sys.argv[1:])

