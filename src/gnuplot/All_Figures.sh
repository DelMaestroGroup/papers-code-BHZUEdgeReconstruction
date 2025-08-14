gnuplot Combined_Akyw_plot.gp
pdflatex OS_resolved_Akyw_for_Vo_9p5_36x36_SSW.tex
pdftoppm -r 700 -png OS_resolved_Akyw_for_Vo_9p5_36x36_SSW.pdf OS_resolved_Akyw_for_Vo_9p5_36x36_SSW
mv OS_resolved_Akyw_for_Vo_9p5_36x36_SSW-1.png ../../figures/OS_resolved_Akyw_for_Vo_9p5_36x36_SSW.png

gnuplot Plot_N_vs_V0.gp
pdflatex Tot_N_vs_V0_for_12x3_GCE.tex
pdftoppm -r 700 -png Tot_N_vs_V0_for_12x3_GCE.pdf Tot_N_vs_V0_for_12x3_GCE
mv Tot_N_vs_V0_for_12x3_GCE-1.png Tot_N_vs_V0_for_12x3_GCE.png

gnuplot TriplePlot_GSE_vs_Sz.gp
pdflatex GSE_vs_Sz_for_diff_Vo_12x3_triple.tex
pdftoppm -r 700 -png GSE_vs_Sz_for_diff_Vo_12x3_triple.pdf GSE_vs_Sz_for_diff_Vo_12x3_triple
mv GSE_vs_Sz_for_diff_Vo_12x3_triple-1.png GSE_vs_Sz_for_diff_Vo_12x3_triple.png

gnuplot Plot_del_ni.gp
pdflatex del_ni_12x3_w4.tex
pdftoppm -r 700 -png del_ni_12x3_w4.pdf del_ni_12x3_w4
mv del_ni_12x3_w4-1.png del_ni_12x3_w4.png

gnuplot Plot_szi.gp
pdflatex szi_12x3_w4.tex
pdftoppm -r 700 -png szi_12x3_w4.pdf szi_12x3_w4
mv szi_12x3_w4-1.png szi_12x3_w4.png

gnuplot Comparison_Akw_plot.gp
pdflatex Akyw_for_diff_eps_Vo_10p0_36x36_SSW.tex
pdftoppm -r 700 -png Akyw_for_diff_eps_Vo_10p0_36x36_SSW.pdf Akyw_for_diff_eps_Vo_10p0_36x36_SSW
mv Akyw_for_diff_eps_Vo_10p0_36x36_SSW-1.png Akyw_for_diff_eps_Vo_10p0_36x36_SSW.png

gnuplot Plot_avg_nrx.gp
pdflatex Avg_nrx_vs_rx_for_12x3_DMRG.tex
pdftoppm -r 700 -png Avg_nrx_vs_rx_for_12x3_DMRG.pdf Avg_nrx_vs_rx_for_12x3_DMRG
mv Avg_nrx_vs_rx_for_12x3_DMRG-1.png Avg_nrx_vs_rx_for_12x3_DMRG.png

gnuplot Plot_FS_exchange_vs_1byLx.gp
pdflatex FS_exchange_vs_1byLx_Nx3_DMRG.tex
pdftoppm -r 700 -png FS_exchange_vs_1byLx_Nx3_DMRG.pdf FS_exchange_vs_1byLx_Nx3_DMRG
mv FS_exchange_vs_1byLx_Nx3_DMRG-1.png FS_exchange_vs_1byLx_Nx3_DMRG.png

rm -f *.tex *.eps *.aux *.log *-eps-converted-to.pdf
rm Akyw_for_diff_eps_Vo_10p0_36x36_SSW.pdf
rm OS_resolved_Akyw_for_Vo_9p5_36x36_SSW.pdf

mv *.pdf ../../figures/
mv *.png ../../figures/

cd ../lattice_plot
pdflatex Confining_Lattice_new.tex
pdftoppm -r 700 -png Confining_Lattice_new.pdf Confining_Lattice_new
mv Confining_Lattice_new-1.png Confining_Lattice_new.png

rm *.log *.aux 
mv *.pdf ../../figures/
mv *.png ../../figures/
