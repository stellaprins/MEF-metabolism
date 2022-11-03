Step2 notes:


mefAtpASEAnalysis takes 'dataRaw' and 'dataKcat' files from step 1 and generates derived datasets capturing different classes of proteins according to EnerSysGo classification via the creation of a data summary array

Use same compression as used in FBA but don't have to run FBA itself. The 'squeeze' or compression vector has to be the same though.

The script then takes as a specific example the tRNA-AA ligases and plots the result.

Finally, summary tables of different classes of protein abundances that are compressed a classified according to whether or not they use ATP etc are stored.



DATA
'dataRaw'
'dataKcat'

'Process_reactions_working2.xlsx'

FUNCTIONS
mapExpression2Matrix

