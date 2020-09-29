set title 'Theory compared to 20 billion random alignments'
set xlabel 'Predicted Probability'
set ylabel 'Observed Probability'
set format "10^{%L}"
set log
set zero 0
set tics out

#set term postscript landscape 18	# for sideways, full page
set term pdf size 3,3
set output "gnuplot.pdf"
#set term dumb 125 37
#set term postscript portrait 22	# for inclusion in LaTeX documents.
#set size 1,1		# desired aspect ratio

#set term gif size 1024,768
#set output "gnuplot.gif"
plot  "AT-HS.NOSEQ.eval" using 6:7 notitle with d, "AT-HS.allGO.eval" using 6:7 notitle with d, "CE-AT.NOSEQ.eval" using 6:7 notitle with d, "CE-AT.allGO.eval" using 6:7 notitle with d, "CE-DM.NOSEQ.eval" using 6:7 notitle with d, "CE-DM.allGO.eval" using 6:7 notitle with d, "CE-HS.NOSEQ.eval" using 6:7 notitle with d, "CE-HS.allGO.eval" using 6:7 notitle with d, "CE-MM.NOSEQ.eval" using 6:7 notitle with d, "CE-MM.allGO.eval" using 6:7 notitle with d, "CE-SC.NOSEQ.eval" using 6:7 notitle with d, "CE-SC.allGO.eval" using 6:7 notitle with d, "DM-AT.NOSEQ.eval" using 6:7 notitle with d, "DM-AT.allGO.eval" using 6:7 notitle with d, "DM-HS.NOSEQ.eval" using 6:7 notitle with d, "DM-HS.allGO.eval" using 6:7 notitle with d, "MM-AT.NOSEQ.eval" using 6:7 notitle with d, "MM-AT.allGO.eval" using 6:7 notitle with d, "MM-DM.NOSEQ.eval" using 6:7 notitle with d, "MM-DM.allGO.eval" using 6:7 notitle with d, "MM-HS.NOSEQ.eval" using 6:7 notitle with d, "MM-HS.allGO.eval" using 6:7 notitle with d, "RN-AT.NOSEQ.eval" using 6:7 notitle with d, "RN-AT.allGO.eval" using 6:7 notitle with d, "RN-CE.NOSEQ.eval" using 6:7 notitle with d, "RN-CE.allGO.eval" using 6:7 notitle with d, "RN-DM.NOSEQ.eval" using 6:7 notitle with d, "RN-DM.allGO.eval" using 6:7 notitle with d, "RN-HS.NOSEQ.eval" using 6:7 notitle with d, "RN-HS.allGO.eval" using 6:7 notitle with d, "RN-MM.NOSEQ.eval" using 6:7 notitle with d, "RN-MM.allGO.eval" using 6:7 notitle with d, "RN-SC.NOSEQ.eval" using 6:7 notitle with d, "RN-SC.allGO.eval" using 6:7 notitle with d, "RN-SP.NOSEQ.eval" using 6:7 notitle with d, "RN-SP.allGO.eval" using 6:7 notitle with d, "SC-AT.NOSEQ.eval" using 6:7 notitle with d, "SC-AT.allGO.eval" using 6:7 notitle with d, "SC-DM.NOSEQ.eval" using 6:7 notitle with d, "SC-DM.allGO.eval" using 6:7 notitle with d, "SC-HS.NOSEQ.eval" using 6:7 notitle with d, "SC-HS.allGO.eval" using 6:7 notitle with d, "SC-MM.NOSEQ.eval" using 6:7 notitle with d, "SC-MM.allGO.eval" using 6:7 notitle with d, "SP-AT.NOSEQ.eval" using 6:7 notitle with d, "SP-AT.allGO.eval" using 6:7 notitle with d, "SP-CE.NOSEQ.eval" using 6:7 notitle with d, "SP-CE.allGO.eval" using 6:7 notitle with d, "SP-DM.NOSEQ.eval" using 6:7 notitle with d, "SP-DM.allGO.eval" using 6:7 notitle with d, "SP-HS.NOSEQ.eval" using 6:7 notitle with d, "SP-HS.allGO.eval" using 6:7 notitle with d, "SP-MM.NOSEQ.eval" using 6:7 notitle with d, "SP-MM.allGO.eval" using 6:7 notitle with d, "SP-SC.NOSEQ.eval" using 6:7 notitle with d, "SP-SC.allGO.eval" using 6:7 notitle with d
