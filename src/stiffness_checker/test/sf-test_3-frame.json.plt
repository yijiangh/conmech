# CONMECH FRAME STRUCTUAL ANALYSIS RESULTS# sf-test_3-frame.json
# Sun Aug  7 05:07:24 2112
# G N U P L O T   S C R I P T   F I L E 
set autoscale
unset border
set pointsize 1.0
set xtics; set ytics; set ztics; 
unset zeroaxis
unset key
unset label
set size ratio -1    # 1:1 2D axis scaling 
# set view equal xyz # 1:1 3D axis scaling 
# NODE NUMBER LABELS
set label ' 1' at   7.0000e-01,  -6.5052e-11,   0.0000e+00
set label ' 2' at  -3.5000e-01,   6.0622e-01,   0.0000e+00
set label ' 3' at  -3.5000e-01,  -6.0622e-01,   0.0000e+00
set label ' 4' at  -2.1259e-09,  -6.5052e-11,   1.0000e+00
# ELEMENT NUMBER LABELS
set label ' 1' at   3.5000e-01,  -6.5052e-11,   5.0000e-01
set label ' 2' at  -1.7500e-01,   3.0311e-01,   5.0000e-01
set label ' 3' at  -1.7500e-01,  -3.0311e-01,   5.0000e-01
  set parametric
  set view 60, 70,  1.00 
  set view equal xyz # 1:1 3D axis scaling 
  unset key
  set xlabel 'x'
  set ylabel 'y'
  set zlabel 'z'
set title "sf-test_3-frame.json"
unset clip; 
set clip one; set clip two
  splot '/Users/yijiangh/Dropbox (MIT)/Projects/conmech/conmech/src/stiffness_checker/test/sf-test_3-frame.json_original.dat' using 2:3:4 title 'original shape' with linespoints  linewidth 1 linetype 1 pointtype 6 
