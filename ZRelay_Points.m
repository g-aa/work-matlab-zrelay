function [ TF ] = ZRelay_Points( Zset, Pv )
    
    % параметры характеристик:
    Line_num = [3; 4; 8;];  % число пр€мых дл€ каждой из характеристик
    
    for i = 1:1:Zset.NUM

        Lset.Lines{i} = Line_num(i);       % число пр€мы в характеристике
        Lset.phi{i} = zeros(Lset.NUM, 1);  % углы наклона пр€мых
        Lset.R0{i} = zeros(Lset.NUM, 1);   % смешенние по оси R
        Lset.X0{i} = zeros(Lset.NUM, 1);   % смещение по оси X
        Lset.kf{i} = zeros(Lset.NUM, 1);   % коэффициент наклона отсь R
        Lset.Ps{i} = zeros(Lset.NUM, 2);   % число точек

    end
    


% пр€ма€ 1:
Lset.phi(1,:) = (Zset.psi - 5);
Lset.R0(1,:) = Zset.RL;
Lset.X0(1,:) = 0;

% пр€ма€ 2:
Lset.phi(2,:) = - 5;
Lset.R0(2,:) = 0;
Lset.X0(2,:) = Zset.ZL*sind(Zset.psi) + Zset.ZL*cosd(Zset.psi)*tand(Lset.phi(2,:));

% пр€ма€ 3:
Lset.phi(3,:) = 105;
Lset.R0(3,:) = -Zset.RL/8;
Lset.X0(3,:) = 0;

% пр€ма€ 4:
Lset.phi(4,:) = -5;
Lset.R0(4,:) = 0;
Lset.X0(4,:) = (Zset.Zs./100)*sind(180 + Zset.psi);


% расчет точек пересечени€ характеристики:
Lset.kf = tand(Lset.phi);
[ Lset.Ps ] = points(Lset.kf, Lset.R0, Lset.X0);

% ===:
[ ~ ] = ZRelayMD(Lset.Ps, Pv);

% ===:
[ ~ ] = ZRelaySD(Lset.Ps, Pv);


%figure;
hold on;
for i = 1:1:Lset.NUM
    if (i < Lset.NUM)
        line([Lset.Ps(i,1), Lset.Ps(i+1,1)], [Lset.Ps(i,2), Lset.Ps(i+1,2)], 'LineWidth', 2);
    else
        line([Lset.Ps(i,1), Lset.Ps(Lset.NUM-i+1,1)], [Lset.Ps(i,2), Lset.Ps(Lset.NUM-i+1,2)], 'LineWidth', 2);
    end
end

% оси:
line([-1.5*Zset.ZL, 1.5*Zset.ZL], [0, 0], 'Color', 'k', 'LineWidth', 1);
line([0, 0], [-1.5*Zset.ZL, 1.5*Zset.ZL], 'Color', 'k', 'LineWidth', 1);

% радиус вектор:
line([0, Pv(1)], [0, Pv(2)], 'Color', 'r', 'LineWidth', 2, 'Marker', '*');
grid minor;
hold off;
end


% расчет точек  характеристики ===========================================
function [ ps ] = points(kf, r0, x0)
    tic;
    DIM = size(kf, 1);
    ps = zeros(DIM ,2);
        
    for i = 1:1:DIM
        if (i < DIM)
            A = [-kf(i,:), 1; -kf(i+1,:), 1];
            b = [x0(i,:) - kf(i,:)*r0(i,:); 
                 x0(i+1,:) - kf(i+1,:)*r0(i,:)];
        else
            A = [-kf(i,:), 1; -kf(DIM-i+1,:), 1];
            b = [x0(i,:) - kf(i,:)*r0(i,:); 
                 x0(DIM-i+1,:) - kf(DIM-i+1,:)*r0(DIM-i+1,:)];
        end 
        ps(i,:) = A\b;
    end
    disp(['¬рем€ работы points: ',num2str(toc*1000), 'мс']);
end


% расчет момента характеристики ƒ« =======================================



% матрица с поиском детерминантом ========================================
function [ TF ] = ZRelayMD(Ps, Pk)
    tic;
    TF = true;  % реле сработало
    DIM = size(Ps, 1);
    
    for i =1:1:DIM
        
        if (i < DIM)
            M = [Ps(i,:), 1; Ps(i+1,:), 1; Pk, 1];
        else
            M = [Ps(i,:), 1; Ps(DIM-i+1,:), 1; Pk, 1];
        end 
        
        if (det(M) <= 0)
            TF = false; % реле не сработало
            break;
        end    
    end
    disp(['¬рем€ расчета MD: ',num2str(toc*1000), 'мс']);
    disp(['—осто€ние реле: ', num2str(TF)]);
end


% статический детерминант ================================================
function [ TF ] = ZRelaySD(Ps, Pk)
    tic;
    TF = true;  % реле сработало
    DIM = size(Ps, 1);
    for i =1:1:DIM  
        if (i < DIM)
            detM = Ps(i,1)*(Ps(i+1,2)-Pk(2)) + Ps(i+1,1)*(Pk(2)-Ps(i,2)) + Pk(1)*(Ps(i,2)-Ps(i+1,2));
        else
            detM = Ps(i,1)*(Ps(DIM-i+1,2)-Pk(2)) + Ps(DIM-i+1,1)*(Pk(2)-Ps(i,2)) + Pk(1)*(Ps(i,2)-Ps(DIM-i+1,2));
        end 
      
        if (detM <= 0)
            TF = false; % реле не сработало
            break;
        end    
    end
    disp(['¬рем€ расчета SD: ',num2str(toc*1000), 'мс']);
    disp(['—осто€ние реле: ', num2str(TF)]);
end