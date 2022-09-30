# Unit Commitment solved with/without Benders Decomposition
- This is a MATLAB program example for Unit Commitment problem. You could choose to solve it with/without Benders Decomposition (BD). 
- You are suggested to read "[UC.pdf](https://github.com/a280558071/UC_Solvedby_BendersDecomp/blob/main/UC.pdf)" (in Chinese) for more information. (Thanks to Biao Zhao (赵彪) in [my research group](https://xinweishen.com/group.html) for his contribution in writing this document)
- The YALMIP function export() is used to export optimization model parameters for BD. 
# Case data and model formulations
- are given in  "[UC.pdf](https://github.com/a280558071/UC_Solvedby_BendersDecomp/blob/main/UC.pdf)" and [1].
- [1] G. Morales-España, J. M. Latorre and A. Ramos, "Tight and Compact MILP Formulation for the Thermal Unit Commitment Problem," in IEEE Transactions on Power Systems, vol. 28, no. 4, pp. 4897-4908, Nov. 2013, doi: 10.1109/TPWRS.2013.2251373.
# Must-include package
- [**YALMIP**](https://yalmip.github.io/) 
- Both ".m" files are calling [**GUROBI**](https://www.gurobi.com/) to solve the problem, but other solvers (e.g. **CPLEX**) could be used as well. (see what solvers they support in abovementioned link for YALMIP).
 # Note
- Sth. would happen if you try to include/exclude Cons. (6) when solve it with BD——The total iteration numbers of BD would change. Just try it by yourself! 

