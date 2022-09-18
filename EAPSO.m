function [BestCost,BestValue] =EAPSO(fhd,nPop,nVar,VarMin,VarMax,MaxIt,X)
MaxIt=floor(MaxIt/nPop*2);
VarSize=[1 nVar];  
empty_particle.Position=[];
empty_particle.Cost=[];
empty_particle.Velocity=[];
empty_particle.Best.Position=[];
empty_particle.Best.Cost=[];
particle=repmat(empty_particle,nPop,1);
GlobalBest.Cost=inf;
for i=1:nPop
    particle(i).Position=X(i,:);
    particle(i).Velocity=zeros(VarSize);
    particle(i).Cost= fhd(particle(i).Position);
    particle(i).Best.Position=particle(i).Position;
    particle(i).Best.Cost=particle(i).Cost;
    if particle(i).Best.Cost<GlobalBest.Cost
        GlobalBest=particle(i).Best;
    end
end
pcount=2;
for i=1:nPop
    cosp(i)=particle(i).Cost;
end
[~,inp]=sort(cosp);
PART(1)=particle(inp(1)).Best;
PART(2)=particle(inp(2)).Best;
NP=nPop;
BestCost(1)=min(cosp);
GART(1)=GlobalBest;
gcount=1;
VelMax=(VarMax-VarMin);
VelMin=-VelMax;
for it=2:MaxIt
    
    for i=1:nPop
        cosp(i)=particle(i).Best.Cost;
    end
    [~,ind]=sort(cosp);
    WINNER=ind(1:nPop/2);
    LOSER= ind(nPop/2+1:nPop);
    for i=1:nPop/2
        a=randperm(length(PART),1);
        b=randperm(length(PART),1);
        while a==b
            a=randperm(length(PART),1);
            b=randperm(length(PART),1);
        end
        if PART(a).Cost>PART(b).Cost
            a=b;
        end
        if length(GART)>1
            c=randperm(length(GART),1);
            d=randperm(length(GART),1);
            while c==d
                c=randperm(length(GART),1);
                d=randperm(length(GART),1);
            end
            if GART(c).Cost>GART(d).Cost
                b=d;
            else
                b=c;
            end
        else
            b=1;
        end
        c=randperm(nPop/2,1);
        q=randperm(nPop/2,1);
        while c==q
            c=randperm(nPop/2,1);
            q=randperm(nPop/2,1);
        end
        if  particle(WINNER(q)).Best.Cost<particle(WINNER(c)).Best.Cost
            c=q;
        end
        w=rand(1,nVar);
        MM=particle(LOSER(i)).Position;
        F1=rand(1,nVar);
        F2=rand(1,nVar);
        if cosp(LOSER(i))<mean(cosp(LOSER))
            if PART(a).Cost<GART(b).Cost && PART(a).Cost<particle(WINNER(c)).Best.Cost 
                particle(LOSER(i)).Velocity=w.*particle(LOSER(i)).Velocity+F1.*(PART(a).Position-MM)+F2.*(GlobalBest.Position-MM);
            elseif PART(a).Cost>GART(b).Cost && GART(b).Cost<particle(WINNER(c)).Best.Cost
                particle(LOSER(i)).Velocity=w.*particle(LOSER(i)).Velocity+F1.*(GART(b).Position-MM)+F2.*(GlobalBest.Position-MM);
            elseif  PART(a).Cost>particle(WINNER(c)).Best.Cost &&  GART(b).Cost>particle(WINNER(c)).Best.Cost
                particle(LOSER(i)).Velocity=w.*particle(LOSER(i)).Velocity+F1.*(particle(WINNER(c)).Best.Position-MM)+F2.*(GlobalBest.Position-MM);
            end
        else
             if PART(a).Cost>GART(b).Cost && PART(a).Cost>particle(WINNER(c)).Best.Cost 
                particle(LOSER(i)).Velocity=w.*particle(LOSER(i)).Velocity+F1.*(GART(b).Position-MM)+F2.*(particle(WINNER(c)).Best.Position-particle(LOSER(i)).Position);
             elseif PART(a).Cost<GART(b).Cost && GART(b).Cost>particle(WINNER(c)).Best.Cost
                particle(LOSER(i)).Velocity=w.*particle(LOSER(i)).Velocity+F1.*(PART(a).Position-MM)+F2.*(particle(WINNER(c)).Best.Position-particle(LOSER(i)).Position);
             elseif PART(a).Cost<particle(WINNER(c)).Best.Cost && GART(b).Cost<particle(WINNER(c)).Best.Cost
                 particle(LOSER(i)).Velocity=w.*particle(LOSER(i)).Velocity+F1.*(PART(a).Position-MM)+F2.*(GART(b).Position-particle(LOSER(i)).Position); 
             end
        end
        particle(LOSER(i)).Velocity = max(particle(LOSER(i)).Velocity,VelMin);
        particle(LOSER(i)).Velocity = min(particle(LOSER(i)).Velocity,VelMax);
        particle(LOSER(i)).Position = particle(LOSER(i)).Position + particle(LOSER(i)).Velocity;
        IsOutside=(particle(LOSER(i)).Position<VarMin | particle(LOSER(i)).Position>VarMax);
        particle(LOSER(i)).Velocity(IsOutside)=-particle(LOSER(i)).Velocity(IsOutside);
        particle(LOSER(i)).Position = max(particle(LOSER(i)).Position,VarMin);
        particle(LOSER(i)).Position = min(particle(LOSER(i)).Position,VarMax);
        particle(LOSER(i)).Cost = fhd(particle(LOSER(i)).Position);
        if particle(LOSER(i)).Cost<particle(LOSER(i)).Best.Cost
            particle(LOSER(i)).Best.Position=particle(LOSER(i)).Position;
            particle(LOSER(i)).Best.Cost=particle(LOSER(i)).Cost;
            pcount=pcount+1;
            if pcount<=NP
                PART(pcount).Position=particle(LOSER(i)).Position;
                PART(pcount).Cost=particle(LOSER(i)).Cost;
            else
                a=randperm(NP,1);
                b=randperm(NP,1);
                while a==b
                    a=randperm(NP,1);
                    b=randperm(NP,1);
                end
                if PART(a).Cost<PART(b).Cost
                    if particle(LOSER(i)).Cost<PART(b).Cost
                        PART(b).Cost=particle(LOSER(i)).Cost;
                        PART(b).Position=particle(LOSER(i)).Position;
                    end
                else
                    if particle(LOSER(i)).Cost<PART(a).Cost
                        PART(a).Cost=particle(LOSER(i)).Cost;
                        PART(a).Position=particle(LOSER(i)).Position;
                    end
                end
            end
        end
        if particle(LOSER(i)).Best.Cost<GlobalBest.Cost
            GlobalBest=particle(LOSER(i)).Best;
        end
    end
    gcount=gcount+1;
    if gcount<=NP
        GART(gcount).Cost= GlobalBest.Cost;
        GART(gcount).Position= GlobalBest.Position;
    else
        a=randperm(length(GART),1);
        b=randperm(length(GART),1);
        while a==b
            a=randperm(length(GART),1);
            b=randperm(length(GART),1);
        end
        if GART(b).Cost>GART(a).Cost
            a=b;
        end
        if  GlobalBest.Cost<GART(a).Cost
            GART(a).Cost= GlobalBest.Cost;
            GART(a).Position= GlobalBest.Position;
        end   
    end
    BestCost(it)=GlobalBest.Cost;
    BestValue=GlobalBest.Cost;
end
end


