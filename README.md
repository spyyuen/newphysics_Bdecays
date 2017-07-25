# semileptonicdecay
Study semileptonic B->Dlnu decays

Help with more run options:<br />
python semileptonicdecay.py -h <br />

To run:
Plot differential decay rate of B->Dlnu decays and R(D) for all cases:<br />
m_l=0<br />
m_l=m_tau<br />
m_l=m_tau with new physics<br />
python semileptonicdecay.py -p -prd <br />

<div align="center">
        <img width="90%" src="plots/differentialDecayRate.png" alt="Run with -p optionscreen" title="Diff decay rate of B->Dlnu"</img>
        <img height="0" width="8px">
        <img width="90%" src="plots/RD_dS.png" alt="Run with -prd option" title="Ratio of decay rates"></img>
        <img height="0" width="8px">
        <img width="90%" src="plots/RD_tanbeta.png" alt="List screen" title="Ratio of decay rates for tanbeta"></img>
</div>


TODO:<br />
Get global constants via scraping<br />
Add legend to plots<br />
Write README
