clc
clear

id = [];
CI_Table = [];
k = 0;
Diff_all = [];
for DGP = [6 1 2 4]
    for Type = [1 3 2]
        k = k+1;
        filename1 = strcat('simulation_DGP',num2str(DGP),'_Type',num2str(Type));   
        load(filename1,'CI_sim','CI_RR','CI_h','CI_pt','T_pre','M','DGP','Au','delta0','CI_YKHS_pt')
        id = [id; [DGP Type max(Au*delta0)]];

        Table_add  =  [median(CI_pt), 0, median(CI_h), 0, median(CI_RR), 0, median(CI_YKHS_pt), 0, median(CI_sim)];
        length     =  [median(CI_pt(:,2))-median(CI_pt(:,1)), median(CI_h(:, 2))-median(CI_h(:, 1)), ...
                       median(CI_RR(:,2))-median(CI_RR(:,1)), median(CI_YKHS_pt(:,2))- median(CI_YKHS_pt(:,1)), median(CI_sim(:,2))-median(CI_sim(:,1))];
        Diff_auc   = length - length(1); 
        Diff_all   = [Diff_all; [DGP Type  Diff_auc] ];
        Diff       = zeros(1,size(Table_add,2));
        Diff(1:3:size(Table_add,2))   =  Diff_auc;
        CI_Table  = [CI_Table; [DGP Type Table_add]; [DGP Type Diff]];
        
    end
end
%%
clc
clear

bound(:,:,1) = [0    0.12;  % 1
                0.3  0.425; % 2
                0.05 0.2;   % 3
                0.3  0.45];  % 4
bound(:,:,2) = [0    4;     % 5
                9    12;    % 6
                0.5  4.5;    % 7
                9    13.5];    % 8
bound(:,:,4) = [0    0.5;     % 5
                0.9    1.2;    % 6
                0.2    0.7;    % 7
                0.9    1.4];    % 8            
bound(:,:,6) = [0    1;  % 1
                2.4  3.5; % 2
                0  1;   % 3
                2.4  3.5];  % 4  
for DGP = [6 1 2 4]
    for Type = [1 3 2]
  
        load(strcat('simulation_DGP',num2str(DGP),'_Type',num2str(Type)),'CI_sim','CI_RR','CI_h','CI_pt','T_pre','M','DGP','Au','delta0','CI_YKHS_pt')

        x = linspace(bound(Type,1,DGP), bound(Type,2,DGP),100);
        rej_h   = 1 - mean((CI_h(:,1) <= x).*(CI_h(:,2) >= x));
        rej_RR  = 1 - mean((CI_RR(:,1) <= x).*(CI_RR(:,2) >= x));       
        rej_YKHS_pt = 1 - mean((CI_YKHS_pt(:,1) <= x).*(CI_YKHS_pt(:,2) >= x));
        rej_sim = 1 - mean((CI_sim(:,1) <= x).*(CI_sim(:,2) >= x));

        compareRR(DGP,Type) = max(rej_h - rej_RR);
        compareYKHS(DGP,Type) = max(rej_h - rej_YKHS_pt);
        comparesim(DGP,Type) = max(rej_h - rej_sim);
    end
    
end

bound(:,:,1) = [0    0.12;  % 1
                0.3  0.425; % 2
                0.05 0.2;   % 3
                0.3  0.45];  % 4
bound(:,:,2) = [0    4;     % 5
                9    12;    % 6
                0.5  4.5;    % 7
                9    13.5];    % 8
bound(:,:,4) = [0    0.5;     % 5
                0.9    1.2;    % 6
                0.2    0.7;    % 7
                0.9    1.4];    % 8            
bound(:,:,6) = [0    1;  % 1
                2.4  3.5; % 2
                0  1;   % 3
                2.4  3.5];  % 4            
k = 0;
for DGP = [6 1 2 4]
    for Type = [1 3 2]
        k = k+1;
        load(strcat('simulation_DGP',num2str(DGP),'_Type',num2str(Type)),'CI_sim','CI_RR','CI_h','CI_pt','T_pre','M','DGP','Au','delta0','CI_YKHS_pt')
   
        figure(k)

        x = linspace(bound(Type,1,DGP), bound(Type,2,DGP),100);
        
%         subplot(1,2,1)
        rej_RR  = 1 - mean((CI_RR(:,1) <= x).*(CI_RR(:,2) >= x));
        rej_h   = 1 - mean((CI_h(:,1) <= x).*(CI_h(:,2) >= x));
        
        x2 = linspace(bound(Type,1,DGP), bound(Type,2,DGP),25);
        rej_YKHS_pt = 1 - mean((CI_YKHS_pt(:,1) <= x2).*(CI_YKHS_pt(:,2) >= x2));
        rej_sim = 1 - mean((CI_sim(:,1) <= x2).*(CI_sim(:,2) >= x2));

        plot(x, rej_h,  'Color','red','LineWidth',1), hold on
        plot(x, rej_RR, 'Color','blue','LineStyle','--','LineWidth',1.1), hold on
        plot(x2, rej_YKHS_pt,'g-o','LineWidth',1), hold on
        plot(x2, rej_sim,'k-.','LineWidth',1), hold on

        t = (x2(2) - x2(1));
        plot(bound(Type,:,DGP): t: max(Au*delta0),0.1*ones(size(bound(Type,:,DGP): t: max(Au*delta0))),'k*')
        xlim(bound(Type,:,DGP))
        grid on
        if DGP==6 && Type == 1
            legend('Modified C','RR23','YKHS','Simple', '\Theta_I','Location','Northwest')
        end
        
%         subplot(1,2,2)
%         rej_RR  = 1 - mean((CI_RR(:,1) <= -x).*(CI_RR(:,2) >= -x));
%         rej_h   = 1 - mean((CI_h_adj(:,1) <= -x).*(CI_h_adj(:,2) >= -x));
%         
%         rej_YKHS_pt = 1 - mean((CI_YKHS_pt(:,1) <= -x2).*(CI_YKHS_pt(:,2) >= -x2));
%         rej_sim = 1 - mean((CI_sim(:,1) <= -x2).*(CI_sim(:,2) >= -x2));
% 
%         plot( -x, rej_h,  'Color','red','LineWidth',1), hold on
%         plot( -x, rej_RR, 'Color','blue','LineStyle','--','LineWidth',1.1), hold on
%         plot(-x2, rej_YKHS_pt,'g-o','LineWidth',1), hold on
%         plot(-x2, rej_sim,'k-.','LineWidth',1), hold on
% 
%         t = (x2(2) - x2(1));
%         plot(-bound(Type,1,DGP): -t: -max(Au*delta0),0.1*ones(size(bound(Type,:,DGP): t: max(Au*delta0))),'k*')
%         xlim([-bound(Type,2,DGP), -bound(Type,1,DGP)])
%         grid on
        saveas(gcf, strcat('DGP',num2str(DGP),'_Type',num2str(Type)), 'png');
    end
    
end

