%% Coded by
% Mohamed Mohamed El-Sayed Atyya
% mohamed.atyya94@eng-st.cu.edu.eg




% Euler angles
close all; clear all; clc;
%% constants
G=6.67428E-11; % gravitational constant
open=1;
%% inputs
% taylor series parameters
final=30*24*3600;   % final time of soution
step=3600;            % rough time step
RT=[3.844e8,0,0];    % [R2]
VT=[0,1.023,0]*1e3;              % [V2]
M=[5.972*10^24,1000];       % [M1,M2]
mu=G*(sum(M));
alpha=0; 
beta=0;
gamma=0;
alpha_old=10; 
beta_old=10;
gamma_old=10;
%% solution
% solution by four terms of taylor series
tl4=step;
for l4=2:length(0:step:final)
    zero=RT(l4-1,1:3);
    first=VT(l4-1,1:3);
    second=-mu*RT(l4-1,1:3)/(norm(RT(l4-1,1:3)))^3;
    third=0*mu/(norm(RT(l4-1,1:3)))^4*(-norm(RT(l4-1,1:3))*VT(l4-1,1:3)+3*RT(l4-1,1:3)*norm(VT(l4-1,1:3)));
    fourth=-mu*((norm(RT(l4-1,1:3)))^3.*second-3*(norm(RT(l4-1,1:3)))^2*norm(VT(l4-1,1:3))*VT(l4-1,1:3))/(norm(RT(l4-1,1:3)))^6+3*mu*((norm(RT(l4-1,1:3)))^4*(norm(second)*RT(l4-1,1:3)+norm(VT(l4-1,1:3))*VT(l4-1,1:3))-4*(norm(RT(l4-1,1:3)))^3*(norm(VT(l4-1,1:3)))^2*RT(l4-1,1:3))/(norm(RT(l4-1,1:3)))^8;
    RT(l4,1:3)=zero+first*tl4+second*tl4^2/2+third*tl4^3/6+fourth*tl4^4/24;
    VT(l4,1:3)=first+second*tl4+third*tl4^2/2+fourth*tl4^3/6;
end
%% plotting
fig=figure();
menu = uicontrol('Parent',fig,'Style','popupmenu','String',{'121';'212';'313';'131';'232';'323';'123';'231';'312';'123';'213';'321';'Exit'},'Units','centimeters' ,'Position',[17.5,0.25,3,0.5]);
slider_alpha = uicontrol('Parent',fig,'Style','slider','Units','centimeters' ,'Position',[0,0,10,0.5],'value',0,'SliderStep', [0.01,0.1] , 'min',0, 'max',360);
slider_beta = uicontrol('Parent',fig,'Style','slider','Units','centimeters' ,'Position',[0,0.75,10,0.5],'value',0,'SliderStep', [0.01,0.1] , 'min',0, 'max',180);
slider_gamma = uicontrol('Parent',fig,'Style','slider','Units','centimeters' ,'Position',[0,1.5,10,0.5],'value',0,'SliderStep', [0.01,0.1] , 'min',0, 'max',360);
set(gcf,'color','w');
while open==1
    figure(fig);
    cla;
    S = get(menu,'value');
    alpha=get(slider_alpha,'value');
    beta=get(slider_beta,'value');
    gamma=get(slider_gamma,'value');
    if S == 1
        state='121';
    elseif S == 2
        state='212';
    elseif S == 3
        state='313';  % euler
    elseif S == 4
        state='131';
    elseif S == 5
        state='232';
    elseif S == 6
        state='323';
    elseif S == 7
        state='123';  % yaw, pitch and roll
    elseif S == 8
        state='231';
    elseif S == 9
        state='312';
    elseif S == 10
        state='123';
    elseif S == 11
        state='213';
    elseif S == 12
        state='321';
    elseif S == 13
        open=0;
    end
%     if alpha ~= alpha_old || beta ~= beta_old || gamma ~= gamma_old
        [ Q ] = RM ( alpha, beta, gamma, state);
        xx=RT*Q;
%     end
    % Taylor solution
    hold all;
    plot3(RT(:,1),RT(:,2),RT(:,3),'color','b','LineWidth',2);
    plot3(RT(1,1),RT(1,2),RT(1,3),'o','color',[.5,0.2,0.9],'LineWidth',5);
    plot3(xx(:,1),xx(:,2),xx(:,3),'color','r','LineWidth',2);
    plot3(xx(1,1),xx(1,2),xx(1,3),'o','color',[.9,0.2,0.5],'LineWidth',5);
    view(3);
    legend('Original orbit','Original orbit start','Transformed orbit','Transformed orbit start','location','southeastoutside')
    grid on;
    xlabel('X','fontsize',18);
    ylabel('Y','fontsize',18);
    zlabel('Z','fontsize',18);
    title(['\alpha = ' num2str(alpha) '^0, \beta = ' num2str(beta) '^0, \gamma = ' num2str(gamma) '^0'],'fontsize',18);
    alpha_old=alpha; beta_old=beta; gamma_old=gamma;
end
close all;