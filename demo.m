clc
clear all
nPop=100;
nVar=30;
MaxIt=150000;
fhd=@sphere;
VarMin=-100.*ones(1,nVar);
VarMax=100.*ones(1,nVar);
for i=1:nPop
    X(i,:)=VarMin+(VarMax-VarMin).*rand(1,nVar);
end
[BestCost,BestValue] =EAPSO(fhd,nPop,nVar,VarMin,VarMax,MaxIt,X);
plot(1:MaxIt/nPop*2,BestCost,'r')
xlabel('The number of function evaluations','Fontname','Times New Roma','fontsize',15','FontWeight','bold');
ylabel('Fitness value','Fontname','Times New Roma','fontsize',15,'FontWeight','bold');