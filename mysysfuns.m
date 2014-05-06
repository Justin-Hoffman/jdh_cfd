function dvdt = mysysfuns(t,S)

dvdt=zeros(size(S));

%Relevant Variables
s=5*10^-4; %Subsistence salary per capita
p=5*10^-3; %Minimum consumption per capita
Gamma=1*10^-2; %Nature regeneration rate
Lambda=1*10^2; %Nature carrying capacity
Am=1*10^-2; %Normal death rate
AM=7*10^-2; %Famine death rate
Bc=3*10^-2; %Commoner birth rate
Be=3*10^-2; %Elite birth rate

%These two variables, k and Delta, will be varied depending on the scenario
k=0; %Inequality factor. In this case, there are no elites, so k=0
Delta=1.67*10^-5; %Prodution rate

%The following variables are assigned to the rows of array y
Xc=S(1); %Commoner population
Xe=S(2); %Elite population
y=S(3); %Nature "wealth"
w=S(4); %Stockpiled wealth

%(possibly assign following values to y array)
%These variables are required for the Handy system of equations:
Wth=p.*Xc+k*p.*Xe; %Minimum wealth before famine
Cc=min(1,w./Wth)*s.*Xc; %Commoner consumption rate
Ce=min(1,w./Wth)*k*s.*Xe; %Elite consumption rate
Ac=Am+max(0,(1-Cc/(s*Xc)))*(AM-Am); %Commoner death rate
Ae=Am+max(0,(1-Ce/(s*Xe)))*(AM-Am); %Elite death rate

dvdt(1)=Bc*Xc-Ac*Xc; %Rate of change of commoner population
dvdt(2)=Be*Xe-Ae*Xe; %Rate of change of elite population
dvdt(3)=Gamma*y*(Lambda-y)-Delta*Xc*y; %Rate of change of nature's resources
dvdt(4)=Delta*Xc*y-Cc-Ce; %Rate of change of stockpiled wealth

end