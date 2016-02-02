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

def format_desc(desc) :
  lines = []
  ns = ''
  ts = ''
  paras = desc.split("\n\n")
  for para in paras :
    ts = ''
    ls = para.split(" ")
    for i,word in enumerate(ls) :
      ts += word + ' '
      if len(ts) > 70 or i == len(ls) - 1 :
        lines.append(ts[:ts.rfind(" ")])
        ts = ts[ts.rfind(" ") + 1:]
    lines.append("")
  return '\n'.join(lines)

def run(args) :
  # defaults
  pos_thresh = 3
  neg_thresh = -3
  n_points = 4
  # process args
  desc = "This program identifies potential purine decoys in DNA. In this context, a purine decoy is one that may need to switch from anti to cis or vise vesa. A potential decoy is identified by having difference density peaks in expected locations around the purine. We look for difference ensity peaks around 5  expecte coordiniates. The expected locations for negative peaks are near the atoms N1, C2 (both on the WC edge), and C8 (opposite the WC edge). The expected location for positive peaks are near C4 and N7 but out from the base, in the plane thereof.\n\nAfter identifying which of the 5 expected coordiniates have difference density near them (if any), we categorize the potential decoy into 1 of three strength category. Super strong decoy candidates have a negative peak near C8 and 2 additional peaks. Strong decoy candidates are not super strong that have either  positive peaks near BOTH positive expected coordinates (near C4 and N7 but out from the base) OR has at least three of the five ponts with peaks near by. Weak decoys has at least two of the five ponts with peaks near by. "


  parser = argparse.ArgumentParser(description=format_desc(desc),
                                 formatter_class=argparse.RawTextHelpFormatter)
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
  parser.add_argument('-k','--write_kin',help=hs,action='store_true')
  parser.add_argument('--keep_files',help='Keep generated pdb files',
                      action='store_true')
  hs = 'Keep output simple'
  parser.add_argument('--simple_out',help=hs,action='store_true')
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
    print "#"*21,args.out_file_name
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
  ddab.set_potential_decoys(log=log)
  if not args.keep_files: ddab.clean_up_files(except_these=[args.out_file_name])
  if not args.simple_out : ddab.write_potential_decoy_summary(log=log)

  if args.out_file_name : log.close()

if __name__ == '__main__' :
  run(sys.argv[1:])

