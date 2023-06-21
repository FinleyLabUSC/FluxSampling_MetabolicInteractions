costlessGrowth('Results/Aero_All','Models/threeModels.mat','Medium/CSources.mat','Medium/minMed.mat',2,1,1,'all')
costlessGrowth('Results/Anaero_All','Models/threeModels.mat','Medium/CSources.mat','Medium/minMed.mat',2,1,0,'all')
processCostlessData('Results/M_all',{'Results/Aero_All','Results/Anaero_All'},{'aer','anaer'})

costlessGrowth('Results/Aero_Most','Models/threeModels.mat','Medium/CSources.mat','Medium/minMed.mat',2,1,1,'most')
costlessGrowth('Results/Anaero_Most','Models/threeModels.mat','Medium/CSources.mat','Medium/minMed.mat',2,1,0,'most')
processCostlessData('Results/M_Most',{'Results/Aero_Most','Results/Anaero_Most'},{'aer','anaer'})

costlessGrowth('Results/Aero_Any','Models/threeModels.mat','Medium/CSources.mat','Medium/minMed.mat',2,1,1,'any')
costlessGrowth('Results/Anaero_Any','Models/threeModels.mat','Medium/CSources.mat','Medium/minMed.mat',2,1,0,'any')
processCostlessData('Results/M_Any',{'Results/Aero_Any','Results/Anaero_Any'},{'aer','anaer'})