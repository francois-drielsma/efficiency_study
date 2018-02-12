set -e
read -p "Really move to paper?" doit

paper_target="/home/cr67/work/2016-11-18_emittance-analysis/paper/"
paper_source="output/2017-02/"

cp -r ${paper_source}/cuts_summary/*/*.tex ${paper_target}/first_observation_of_emittance_reduction/02-Cuts/
cp -r ${paper_source}/compare_cuts/* ${paper_target}/first_observation_of_emittance_reduction/02-Cuts/Figures
cp -r ${paper_source}/compare_globals ${paper_target}/first_observation_of_emittance_reduction/03-Detectors/Figures

echo "Moved from ${paper_source} to ${paper_target}"