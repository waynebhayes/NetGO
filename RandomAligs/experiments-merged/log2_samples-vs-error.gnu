set title 'Error vs number of samples'
set xlabel 'Number of probability samples'
set ylabel 'Mean fractional error'
set format x "10^{%L}"
set format y "10^{%T}"
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
plot [:][1e-5:2] "log2_samples-vs-error.txt" using 1:3:4 notitle with errorbars
