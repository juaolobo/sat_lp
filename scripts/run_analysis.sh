for i in $(seq 1 10);
do
	echo "running for $i"
	python analysis/analysis_fixing.py uf20-0$i-simplex-feasibility
	python analysis/generate_figs.py uf20-0$i-simplex-feasibility
	python analysis/analysis_fixing.py uf20-0$i-ip-feasibility
	python analysis/generate_figs.py uf20-0$i-ip-feasibility
done
