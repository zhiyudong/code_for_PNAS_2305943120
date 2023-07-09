clear
%D=1; m=0.2; h1=1;

plotFS = true;

%% tunnable parameters
nkx = 200; %60
nky = 200;
nband = 4; % 4 bands in each valley each spin
%D = 50;
%u_a = 0;
%m = gamma1/(2*v^2);

%mu = D - 8;

%target_band = 3; % band where FS lies (count from bottom)
%D_list = linspace(10,200,40); % 10,70,40
u_list = linspace(30,80,100); %30 80,100
n_list = linspace(-2E-3,-7E-3,100); %-2E-3,-7E-3,100 ; % TOTAL carrier density (two spins two valleys)
n_per_isospin_list = n_list/4; % average carrier density PER SPIN PER VALLEY
%BLG h-doped:-1E-3:-1E-4:-6E-3;
%BLG e-doped: 1E-4:1E-4:4E-3;
% 1E-3:2E-4:1E-2; nm^-2
%-1E-2:2E-4:-1E-3;



%% constants
hbar = 6.58*1E-13; %meV s
%aCC = 1.46E-1; % nm
c = 3E8*1E9; %nm/s
me = 0.511*1E9/(c^2);% meV s^2/ nm^2
%m = 0.028*me; % m (meV s^2/ nm^2)
%ma = 0.19*me; % (meV s^2/ nm^2)
%pD = 0.058/aCC;
% v3 = 1.3* 1E5 * 1E9; % nm/s

a = 0.246; % nm
kmax = 0.028*pi/a; %0.03 pi/a is enough for D=100

V = 6000; %meV nm^2
% if TLG, then V=200 for electron doping (n= +1E-3 ~ +1E-2)

% %https://iopscience.iop.org/article/10.1088/0034-4885/76/5/056503/pdf
% gamma0 = 3160;%meV
% gamma1 = 381;
% %gamma2 = -20; %
% gamma3 = 380;
% gamma4 = 140;
% %gamma5 = 38;
% %Delta = 22;
% Deltap = 22;
% deltaAB = 0;

% Parameter below comes from[Jeil Jung and Allan H. MacDonald Phys. Rev. B
% 89, 035405],  i.e. https://journals.aps.org/prb/abstract/10.1103/PhysRevB.89.035405
% also used in https://arxiv.org/pdf/2302.00682.pdf
% NOT from  https://arxiv.org/pdf/2303.00742.pdf
% 
gamma0 = 2610;%3160; %meV
gamma1 = 361; %381;%
gamma3 = 283; %380; %
gamma4 = 138;
Deltap =  15; %0; %
deltaAB = 0;

lambda = 0; % Ising SOC %meV

v = sqrt(3)/2 * (a*gamma0/hbar);
v1 = sqrt(3)/2 * (a*gamma1/hbar);
v3 = sqrt(3)/2 * (a*gamma3/hbar);
v4 = sqrt(3)/2 * (a*gamma4/hbar);

% A/B sublattice space
sigma_0 = eye(2);
sigma_x = [0,1;1,0];
sigma_y = [0,-1i;1i,0];
sigma_z = [1,0;0,-1];

% K/K' valley space
tau_0 = eye(2);
tau_x = [0,1;1,0];
tau_y = [0,-1i;1i,0];
tau_z = [1,0;0,-1];



if max(n_list)>0 && min(n_list)>0 % e-dope
    charge_sign = 1;
    %target_bands = 5:6; % two lowest upper bands in two valleys
    conduction_band_K_up = 3; % the lowest upper band in valley K spin up, count from bottom (ascend)
    conduction_band_K_dn = 3; % the lowest upper band in valley K spin dn, count from bottom (ascend)
    conduction_band_Kp_up = 3; % the lowest upper band in valley K' spin up, count from bottom (ascend)
    conduction_band_Kp_dn = 3; % the lowest upper band in valley K' spin dn, count from bottom (ascend)
    sort_direction = 'ascend'; %This will be used to find n lowest energy states in two upper bands ("two" due to two valleys)
elseif max(n_list)<0 && min(n_list)<0 % hole-dope
    charge_sign = -1;
    %target_bands = 3:4; % two highest lower bands in two valleys
    conduction_band_K_up = 2; % the lowest upper band in valley K spin up, count from bottom (ascend)
    conduction_band_K_dn = 2; % the lowest upper band in valley K spin dn, count from bottom (ascend)
    conduction_band_Kp_up = 2; % the lowest upper band in valley K' spin up, count from bottom (ascend)
    conduction_band_Kp_dn = 2; % the lowest upper band in valley K' spin dn, count from bottom (ascend)
    sort_direction = 'descend'; %This will be used to find n highest energy state in two lower bands ("two" due to two valleys)
else
    fprintf("error")
end
%n_list = charge_sign*n_list;
%n_per_spin_list = charge_sign*n_per_spin_list;
% in calculation, always set n to positive, describing e- or hole density
% will map back to real n in Figure

%%
kx_list = kmax*linspace(-1,1,nkx); dkx = kx_list(2)-kx_list(1);
ky_list = kmax*linspace(-1,1,nky); dky = ky_list(2)-ky_list(1);
[kx_mesh,ky_mesh] = meshgrid(kx_list,ky_list);
num_k = nkx*nky;

p = zeros(nkx,nky,2);
epsilonK = zeros(nkx,nky,2);
epsilonKp = zeros(nkx,nky,2);
epsilon = zeros(nkx,nky,nband);
%H = zeros(nband, nband);

%epsilon_VP = zeros(nkx,nky,nband,length(n_list));
%epsilon_IVC = zeros(nkx,nky,nband,length(n_list));
%n_ = zeros(1,nband);
E_kinetic = zeros(1,nband);

epsilon0 = zeros(nkx,nky,nband);
epsilon0K_up = zeros(nkx,nky,nband);
epsilon0K_dn = zeros(nkx,nky,nband);
epsilon0Kp_up = zeros(nkx,nky,nband);
epsilon0Kp_dn = zeros(nkx,nky,nband);
% epsilon_UP = zeros(nkx,nky,nband);
% epsilon_SP = zeros(nkx,nky,nband);
% epsilon_SP_VP = zeros(nkx,nky,nband);
% %epsilon_VP_minus = zeros(nkx,nky,nband);
% epsilon_SP_IVC = zeros(nkx,nky,nband);
% %epsilon_IVC_minus = zeros(nkx,nky,nband);

polarization = zeros(length(n_list),length(u_list));
is_FS_overlap_K_up_and_Kp_dn = zeros(length(n_list),length(u_list));
is_FS_overlap_K_dn_and_Kp_up = zeros(length(n_list),length(u_list));
is_three_pocket_K_up = ones(length(n_list),length(u_list));
is_three_pocket_K_dn = ones(length(n_list),length(u_list));
eh_pockets_num_K_up = ones(length(n_list),length(u_list));
eh_pockets_num_K_dn = ones(length(n_list),length(u_list));
tic
for i_u = 1:length(u_list)
    u = u_list(i_u);
    for i_n = 1:length(n_list)
        fprintf("i_u = %d,\t i_n = %d, \t t=%.2f s\n",i_u,i_n,toc);
        tic
        n = n_list(i_n);
        n_per_isospin = n_per_isospin_list(i_n);
        %         delta_n_list = linspace(0,n,10+1);
        %         delta_n_list(1) = [];
        %         E_tot_SP_VP_list=zeros(1, length(delta_n_list));
        %         E_tot_SP_IVC_list=zeros(1, length(delta_n_list));



        %for i_delta_n = 1:length(delta_n_list)
        %    delta_n = delta_n_list(i_delta_n);
        %fprintf("delta_n=%f\n",delta_n)
        %n_plus  = 1/2 * (n+delta_n);
        %n_minus = 1/2 * (n-delta_n);
        %N_plus = floor(n_plus/(dkx*dky/(2*pi)^2));
        %N_minus = floor(n_minus/(dkx*dky/(2*pi)^2));
        N_total = floor(abs(n)/(dkx*dky/(2*pi)^2));
        %N_each_band = nkx*nky;

        for idx_kx = 1: nkx
            for idx_ky = 1: nky

                %                     if abs(delta_n) > abs(n)
                %                         fprintf("ERROR: delta_n=%f\t, n=%f\n, ",delta_n,n)
                % epsilon0(idx_kx,idx_ky,1:end) = nan(1,nband);
                %                         epsilon_UP(idx_kx,idx_ky,1:end) = nan(1,nband);
                %                         epsilon_SP(idx_kx,idx_ky,1:end) = nan(1,nband);
                %                         epsilon_SP_VP(idx_kx,idx_ky,1:end) = nan(1,nband);
                %                         epsilon_SP_IVC(idx_kx,idx_ky,1:end) = nan(1,nband);
                %                         continue
                %                     end

                px = hbar*kx_list(idx_kx);
                py = hbar*ky_list(idx_ky);
                p(idx_kx,idx_ky,:) = [px,py];

                piK = (px + 1i*py);     piK_dagger = (px - 1i*py);
                piKp = (-px + 1i*py);   piKp_dagger = (-px - 1i*py);


                % Hamiltonian below using [Edward McCann, Mikito Koshino 2013 Rep. Prog. Phys. 76 056503]
                % This hamiltonian differs from that in Eq.26 in [Jeil Jung and Allan H. MacDonald Phys. Rev. B
                % 89, 035405].
                % only on the sign in front of v3
                % In [Jeil...] it is -v3, In [Edwar...] it is v3
                % See Table III in [Jeil...] and discussion below Eq.26

                HK_up =[-u/2+lambda,   v*piK_dagger,           -v4*piK_dagger, v3*piK;
                    v*piK,              -u/2+Deltap+lambda,    gamma1,         -v4*piK_dagger;
                    -v4*piK,            gamma1,                 u/2+Deltap,     v*piK_dagger;
                    v3*piK_dagger,       -v4*piK,               v*piK,          u/2   ];

                %                                     HK =[-u/2+lambda,          v*piK_dagger,   -v4*piK_dagger, v3*piK;
                %                         v*piK,          -u/2+Deltap+lambda,           gamma1,         -v4*piK_dagger;
                %                         -v4*piK,        gamma1,         u/2+Deltap,            v*piK_dagger;
                %                         v3*piK_dagger,  -v4*piK,        v*piK,           u/2   ];

                HK_dn =[-u/2-lambda,    v*piK_dagger,           -v4*piK_dagger, v3*piK;
                    v*piK,              -u/2+Deltap-lambda,     gamma1,         -v4*piK_dagger;
                    -v4*piK,            gamma1,                 u/2+Deltap,     v*piK_dagger;
                    v3*piK_dagger,      -v4*piK,                v*piK,          u/2];


                HKp_up =[-u/2-lambda,  v*piKp_dagger,          -v4*piKp_dagger, v3*piKp;
                    v*piKp,         -u/2+Deltap-lambda,         gamma1,         -v4*piKp_dagger;
                    -v4*piKp,       gamma1,                     u/2+Deltap,     v*piKp_dagger;
                    v3*piKp_dagger, -v4*piKp,                   v*piKp,         u/2   ];

                HKp_dn =[-u/2+lambda,  v*piKp_dagger,          -v4*piKp_dagger, v3*piKp;
                    v*piKp,         -u/2+Deltap+lambda,         gamma1,         -v4*piKp_dagger;
                    -v4*piKp,       gamma1,                     u/2+Deltap,     v*piKp_dagger;
                    v3*piKp_dagger, -v4*piKp,                   v*piKp,         u/2   ];


                %                     H_0    = -1/(2*m) * ((px^2-py^2)*sigma_x + 2*px*py*sigma_y); % for valley K
                %                     H_0_Kp = -1/(2*m) * ((px^2-py^2)*sigma_x - 2*px*py*sigma_y); % for valley K'
                %                     H_w    = v3 * ( px*sigma_x - py*sigma_y) - v3 * a/(4*sqrt(3)*hbar) *  ((px^2-py^2)*sigma_x + 2*px*py*sigma_y);
                %                     H_w_Kp = v3 * (-px*sigma_x - py*sigma_y) - v3 * a/(4*sqrt(3)*hbar) *  ((px^2-py^2)*sigma_x - 2*px*py*sigma_y);
                %                     H_4    = 2*v0*v4/gamma1 * (px^2+py^2)* sigma_0;
                %                     H_4_Kp = 2*v0*v4/gamma1 * (px^2+py^2)* sigma_0;
                %                     H_Delta    = Deltap*v0^2/gamma1^2 * (px^2+py^2)* sigma_0;
                %                     H_Delta_Kp = Deltap*v0^2/gamma1^2 * (px^2+py^2)* sigma_0;
                %                     H_U    = -D * (sigma_z - 2*v0^2/gamma1^2*(px^2+py^2) * sigma_z); % our D = U/2 in that paper
                %                     H_U_Kp = -D * (sigma_z - 2*v0^2/gamma1^2*(px^2+py^2) * sigma_z);
                %                     H_deltaAB    = deltaAB/2 * sigma_z ;
                %                     H_deltaAB_Kp = deltaAB/2 * sigma_z ;

                %H_ch = v0^3 / gamma1^2 * ((px^3 - 3*px*py^2) * sigma_x + (3*px^2*py - py^3) * sigma_y) ;
                %H_s = (delta - 2*v0*v4* (px^2+py^2) / gamma1 ) *sigma_0;
                %H_tr = (gamma2/2 - 2*v0*v3* (px^2+py^2) / gamma1 ) *sigma_x;
                %H_gap = D*(1 - v0^2* (px^2+py^2) / gamma1^2 ) *sigma_z;
                %H_sp = zeros(2,2);%u_a*(1 - 3*v0^2* (px^2+py^2) / gamma1^2 ) *sigma_0;
                %                     HK  = H_0    + H_w    + H_4    + H_Delta    + H_U    + H_deltaAB;
                %                     HKp = H_0_Kp + H_w_Kp + H_4_Kp + H_Delta_Kp + H_U_Kp + H_deltaAB_Kp;

                %                     H0= zeros(16,16);
                %                     H0(1:4,1:4) = HK_up;
                %                     H0(5:8,5:8) = HK_dn;
                %                     H0(9:12,9:12) = HKp_up;
                %                     H0(13:16,13:16) = HKp_dn;
                %
                %                     eval = eigs(H0,size(H0,1));%,nband,'smallestreal'
                %                     eval = sort(eval,'ascend');
                %                     epsilon0(idx_kx,idx_ky,:) = eval;

                eval = eigs(HK_up,size(HK_up,1));
                eval = sort(eval,'ascend');
                epsilon0K_up(idx_kx,idx_ky,:) = eval;

                eval = eigs(HK_dn,size(HK_dn,1));
                eval = sort(eval,'ascend');
                epsilon0K_dn(idx_kx,idx_ky,:) = eval;

                eval = eigs(HKp_up,size(HKp_up,1));
                eval = sort(eval,'ascend');
                epsilon0Kp_up(idx_kx,idx_ky,:) = eval;

                eval = eigs(HKp_dn,size(HKp_dn,1));
                eval = sort(eval,'ascend');
                epsilon0Kp_dn(idx_kx,idx_ky,:) = eval;
            end
        end

        %% calculate chemical potential
        All_states_energy = nan(4*nkx*nky,1);
        A = epsilon0K_up(:,:,conduction_band_K_up);
        All_states_energy(1:nkx*nky , 1) = A(:);
        A = epsilon0K_dn(:,:,conduction_band_K_dn);
        All_states_energy((1*nkx*nky+1):(2*nkx*nky),1 ) = A(:);
        A = epsilon0Kp_up(:,:,conduction_band_Kp_up);
        All_states_energy((2*nkx*nky+1):(3*nkx*nky),1 ) = A(:);
        A = epsilon0Kp_dn(:,:,conduction_band_Kp_dn);
        All_states_energy((3*nkx*nky+1):(4*nkx*nky),1 ) = A(:);

        epsilon0_sort = sort(All_states_energy,sort_direction);
        mu= epsilon0_sort(N_total);
        mu_plane = zeros(size(kx_mesh))+mu;

        %% aaa
        el_filling_K_up = epsilon0K_up(:,:,conduction_band_K_up)<mu; % electron occupation in the lowest upper band of valley K spin up
        el_filling_K_dn = epsilon0K_dn(:,:,conduction_band_K_dn)<mu; % electron occupation in the lowest upper band of valley K spin down
        el_filling_Kp_up = epsilon0Kp_up(:,:,conduction_band_Kp_up)<mu; % electron occupation in the lowest upper band of valley K' spin up
        el_filling_Kp_dn = epsilon0Kp_dn(:,:,conduction_band_Kp_dn)<mu; % electron occupation in the lowest upper band of valley K' spin down
        hole_filling_K_up = ones(nkx,nky) - el_filling_K_up; % hole occupation in the highest lower band of valley K spin up
        hole_filling_K_dn = ones(nkx,nky) - el_filling_K_dn; % hole occupation in the highest lower band of valley K spin dn
        hole_filling_Kp_up = ones(nkx,nky) - el_filling_Kp_up; % hole occupation in the highest lower band of valley K' spin up
        hole_filling_Kp_dn = ones(nkx,nky) - el_filling_Kp_dn; % hole occupation in the highest lower band of valley K' spin dn
        if (charge_sign==1) % e-dope
            FS_overlap_K_up_and_Kp_dn = el_filling_K_up.*el_filling_Kp_dn; %the overlap of between two electron Fermi seas in K and K' (in the higher band)
            FS_overlap_K_dn_and_Kp_up = el_filling_K_dn.*el_filling_Kp_up; %the overlap of between two electron Fermi seas in K and K' (in the higher band)
        elseif (charge_sign==-1) % hole-dope
            FS_overlap_K_up_and_Kp_dn = hole_filling_K_up.*hole_filling_Kp_dn; %the overlap of between two hole Fermi seas in K and K' (in the lower band)
            FS_overlap_K_dn_and_Kp_up = hole_filling_K_dn.*hole_filling_Kp_up; %the overlap of between two hole Fermi seas in K and K' (in the lower band)
        end
        is_FS_overlap_K_up_and_Kp_dn(i_n,i_u) = sum(sum(FS_overlap_K_up_and_Kp_dn))>0;
        is_FS_overlap_K_dn_and_Kp_up(i_n,i_u) = sum(sum(FS_overlap_K_dn_and_Kp_up))>0;

        Count_hole_pockets_K_up = bwconncomp(hole_filling_K_up,4);
        Count_hole_pockets_K_dn = bwconncomp(hole_filling_K_dn,4);
        Count_el_pockets_K_up = bwconncomp(el_filling_K_up,4);
        Count_el_pockets_K_dn = bwconncomp(el_filling_K_dn,4);
        hole_pocket_num_K_up = Count_hole_pockets_K_up.NumObjects;
        hole_pocket_num_K_dn = Count_hole_pockets_K_dn.NumObjects;
        el_pocket_num_K_up = Count_el_pockets_K_up.NumObjects;
        el_pocket_num_K_dn = Count_el_pockets_K_dn.NumObjects;
%         Count_hole_pockets_Kp_up = bwconncomp(hole_filling_Kp_up,4);
%         Count_hole_pockets_Kp_dn = bwconncomp(hole_filling_Kp_dn,4);
%         hole_pocket_num_Kp_up = Count_hole_pockets_Kp_up.NumObjects;
%         hole_pocket_num_Kp_dn = Count_hole_pockets_Kp_dn.NumObjects;
        
        is_three_pocket_K_up(i_n,i_u) = (hole_pocket_num_K_up==3);
        is_three_pocket_K_dn(i_n,i_u) = (hole_pocket_num_K_dn==3);
        eh_pockets_num_K_up(i_n,i_u) = hole_pocket_num_K_up + el_pocket_num_K_up;
        eh_pockets_num_K_dn(i_n,i_u) = hole_pocket_num_K_dn + el_pocket_num_K_dn;


        %% plot band structure
if plotFS==true
    figure(1)
        %%target_bands = 3:4;
%                     for i_band = target_bands
%                         surf(kx_mesh,ky_mesh,epsilon0(:,:,i_band)','edgecolor','none','FaceAlpha',0.5)
%                         hold on
%                     end
% 
%         % plot band dispersion
%                     surf(kx_mesh,ky_mesh,epsilon0K(:,:,conduction_band_K)','edgecolor','none','Facecolor','r','FaceAlpha',0.5)
%                     hold on
%                     surf(kx_mesh,ky_mesh,epsilon0Kp(:,:,conduction_band_Kp)','edgecolor','none','Facecolor','b','FaceAlpha',0.5)
%                     hold on
%                     surf(kx_mesh,ky_mesh,mu_plane,'edgecolor','none','Facecolor',[0.5,0.5,0.5],'FaceAlpha',0.5)
%                     hold off
% 
%         band_top = max(max(epsilon0K(:,:,conduction_band_K)));
%         contour_range = linspace(band_top,band_top-5,10);


%% plot FS

%         subplot(1,2,1)
        contour_range=[mu,mu];
        contour3(kx_mesh,ky_mesh,epsilon0K_up(:,:,conduction_band_K_up)',contour_range,'r','linewidth',2)
        hold on
        contour3(kx_mesh,ky_mesh,epsilon0Kp_dn(:,:,conduction_band_Kp_up)',contour_range,'b','linewidth',2)
        grid off
        hold off
        axis equal
       %xticks([-0.2,0,0.2])
        %yticks([-0.2,0,0.2])
        view(2)
        ax = gca;
        ax.FontSize = 20;

        %xlabel('$k_x [\rm{nm}^{-1}]$','interpreter','latex')
        %ylabel('$k_y [\rm{nm}^{-1}]$','interpreter','latex')
        %zlabel('$\epsilon(k)$','interpreter','latex')
        axis off
        drawnow
% 
%         subplot(1,2,2)
%         contour_range=[mu,mu];
%         contour3(kx_mesh,ky_mesh,epsilon0K_dn(:,:,conduction_band_K_dn)',contour_range,'color',[0.9,0.5,0.1])
%         hold on
%         contour3(kx_mesh,ky_mesh,epsilon0Kp_up(:,:,conduction_band_Kp_dn)',contour_range,'color',[0.3,0.1,0.5])
%         grid off
%         hold off
%         axis equal
% 
%         %%
%         xticks([-0.2,0,0.2])
%         yticks([-0.2,0,0.2])
%         view(2)
%         ax = gca;
%         ax.FontSize = 20;
%         axis off
%         set(gcf, 'Position',  [100, 100, 500, 500])
% 
% 
%         %shading interp
%         xlabel('$k_x [\rm{nm}^{-1}]$','interpreter','latex')
%         ylabel('$k_y [\rm{nm}^{-1}]$','interpreter','latex')
%         zlabel('$\epsilon(k)$','interpreter','latex')
%         drawnow
end

        %% calculate the band structure in the presence of orders

        %                     Psigmaz = 1/2*(sigma_0-charge_sign*sigma_z);
        %                     Ptaux_p = 1/2*(tau_0+tau_x);
        %                     Ptaux_m = 1/2*(tau_0-tau_x);
        %                     Ptauz_p = 1/2*(tau_0+tau_z);
        %                     Ptauz_m = 1/2*(tau_0-tau_z);
        %
        %                     H_UP = H0 - V*(n/4)*kron(tau_0,Psigmaz);
        %
        %                     H_SP = H0 - V*(n/2)*kron(tau_0,Psigmaz);
        %
        %                     H_IVC = H0 - V*(n/2 + delta_n/2)*kron(Ptaux_p,Psigmaz)...
        %                                - V*(n/2 - delta_n/2)*kron(Ptaux_m,Psigmaz);
        %                     %H_IVC = H0 - V*delta_n/2*kron(tau_x,Pz);
        %                     %H_IVC_minus = H0 + V*delta_n/2*kron(tau_x,Pz);
        %
        %                     H_VP = H0 - V*(n/2 + delta_n/2)*kron(Ptauz_p,Psigmaz)...
        %                               - V*(n/2 - delta_n/2)*kron(Ptauz_m,Psigmaz);
        %                     %H_VP = H0 - V*delta_n/2*kron(tau_z,Pz);
        %                     %H_VP_minus = H0 + V*delta_n/2*kron(tau_z,Pz);
        %
        %                     %         eval = eigs(H,nband,'smallestabs');
        %                     %         eval(eval<0) = inf;
        %                     %         eval = sort(eval,'ascend');
        %                     %         epsilon(idx_kx,idx_ky,:) = eval;
        %
        %                     %         [~,eval] = eigs(HK,2,'smallestreal');
        %                     %         [eval,sort_ind] = sort(diag(eval),'ascend');
        %                     %                 epsilonK(idx_kx,idx_ky,:) = eval;
        %                     %
        %                     %                 [~,eval] = eigs(HKp,2,'smallestreal');
        %                     %                 [eval,sort_ind] = sort(diag(eval),'ascend');
        %                     %                 epsilonKp(idx_kx,idx_ky,:) = eval;
        %
        %                     if (i_delta_n==1)
        %                     [~,eval] = eigs(H0,nband,'smallestreal');
        %                     epsilon0(idx_kx,idx_ky,:) = diag(eval);
        %
        %                     [~,eval] = eigs(H_UP,nband,'smallestreal');
        %                     epsilon_UP(idx_kx,idx_ky,:) = diag(eval);
        %
        %                     [~,eval] = eigs(H_SP,nband,'smallestreal');
        %                     epsilon_SP(idx_kx,idx_ky,:) = diag(eval);
        %                     end
        %
        %                     [~,eval] = eigs(H_VP,nband,'smallestreal');
        %                     epsilon_SP_VP(idx_kx,idx_ky,:) = diag(eval);
        %                     %[~,eval] = eigs(H_VP_minus,nband,'smallestreal');
        %                     %epsilon_VP_minus(idx_kx,idx_ky,:) = diag(eval);
        %
        %                     [~,eval] = eigs(H_IVC,nband,'smallestreal');
        %                     epsilon_SP_IVC(idx_kx,idx_ky,:) = diag(eval);
        %                     %[~,eval] = eigs(H_IVC_minus,nband,'smallestreal');
        %                     %epsilon_IVC_minus(idx_kx,idx_ky,:) = diag(eval);
        %                     %[evecs,eval] = eigs(H,nband,'smallestabs');
        %                     %eval(eval<0) = inf;
        %                     % evecs = evecs(:,sort_ind);
        %                 end
        %             end
        %
        %
        %
        %             A = epsilon_UP(:,:,target_bands);
        %             epsilon_UP_sort = sort(A(:),sort_direction);
        %             E_tot_UP = 2*charge_sign* sum(epsilon_UP_sort(1:floor(N/2)))*(dkx*dky/(2*pi)^2)...
        %                 + 4*V/2*(1/16*n^2);
        %
        %             A = epsilon_SP(:,:,target_bands);
        %             epsilon_SP_sort = sort(A(:),sort_direction);
        %             E_tot_SP = charge_sign* sum(epsilon_SP_sort(1:N))*(dkx*dky/(2*pi)^2)...
        %                 + 2*V/2*(1/4*n^2);
        %
        %
        %             A = epsilon_SP_VP(:,:,target_bands);
        %             epsilon_SP_VP_sort = sort(A(:),sort_direction);
        %             %A = epsilon_VP_minus(:,:,4);
        %             %epsilon_VP_minus_sort = sort(A(:));
        %             E_tot_SP_VP_list(i_delta_n) = ...
        %                 charge_sign* sum(epsilon_SP_VP_sort(1:N ))*(dkx*dky/(2*pi)^2) ...
        %                 + 2*V/2*(1/4*n^2 + 1/4*delta_n^2);
        %             %+sum(epsilon_VP_minus_sort(1:N_minus))*(dkx*dky/(2*pi)^2) ...
        %
        %
        %             A = epsilon_SP_IVC(:,:,target_bands);
        %             epsilon_SP_IVC_sort = sort(A(:),sort_direction);
        %             %A = epsilon_IVC_minus(:,:,4);
        %             %epsilon_IVC_minus_sort = sort(A(:));
        %             E_tot_SP_IVC_list(i_delta_n) = ...
        %                 charge_sign* sum(epsilon_SP_IVC_sort(1:N ))*(dkx*dky/(2*pi)^2) ...
        %                 + 2*V/2*(1/4*n^2 + 1/4*delta_n^2);
        %             %+sum(epsilon_IVC_minus_sort(1:N_minus))*(dkx*dky/(2*pi)^2) ...
        %
        %
        % %             A = epsilon0(:,:,target_bands);
        % %             epsilon0_sort = sort(A(:),sort_direction);
        % %             mu= epsilon0_sort(N);
        % %             mu_plane = zeros(size(kx_mesh))+mu;
        % %             figure(1)
        % %             %target_bands = 3:4;
        % %             for i_band = target_bands
        % %                 surf(kx_mesh,ky_mesh,epsilon0(:,:,i_band)','edgecolor','none','FaceAlpha',0.5)
        % %                 hold on
        % %             end
        % %             surf(kx_mesh,ky_mesh,mu_plane,'edgecolor','none','Facecolor',[0.5,0.5,0.5])
        % %             hold off
        % %             %shading interp
        % %             xlabel('$k_x [\rm{nm}^{-1}]$','interpreter','latex')
        % %             ylabel('$k_y [\rm{nm}^{-1}]$','interpreter','latex')
        % %             zlabel('$\epsilon(k)$','interpreter','latex')
        % %             drawnow
        %         end
        %
        %         E_tot_SP_VP_min = nanmin(E_tot_SP_VP_list);
        %         E_tot_SP_IVC_min = nanmin(E_tot_SP_IVC_list);
        %         min_id_SP_VP = find(E_tot_SP_VP_list == E_tot_SP_VP_min);
        %         min_id_SP_IVC = find(E_tot_SP_IVC_list == E_tot_SP_IVC_min);
        %
        %         E_tot = [E_tot_UP, E_tot_SP, E_tot_SP_VP_min,E_tot_SP_IVC_min];
        %
        %         if (E_tot_UP == min(E_tot))
        %             polarization(i_n,i_D) = nan;
        %         elseif (E_tot_SP == min(E_tot))
        %             polarization(i_n,i_D) = 0;
        %         elseif (E_tot_SP_VP_min == nanmin(E_tot))
        %             polarization(i_n,i_D) = min_id_SP_VP/length(delta_n_list);
        %         elseif (E_tot_SP_IVC_min == nanmin(E_tot))
        %             polarization(i_n,i_D) = -min_id_SP_IVC/length(delta_n_list);
        %        end
    end
end

%n_list = charge_sign*n_list;
[n_mesh,u_mesh] = meshgrid(n_list,u_list);

%% plot overlap of Fermi seas in K and K'
figure(2)
% sp(1)=subplot(1,2,1);
% surf(n_mesh,u_mesh,is_FS_overlap_K_up_and_Kp_dn(:,:)','edgecolor','none','FaceAlpha',0.5)
% hold off
% %shading interp
% xlabel('$n_e [\rm{nm}^{-2}]$','interpreter','latex')
% ylabel('$u [\rm{meV}]$','interpreter','latex')
% view(2)
% %mycolormap = customcolormap_preset('red-yellow-blue');
% colorbar;%('southoutside');
% colormap(sp(1),customcolormap_preset('red-yellow-blue'));
% caxis_max = 1;
% %clim([-caxis_max caxis_max])
% grid off
% ax = gca;
% ax.FontSize = 20;
% 
% sp(2) = subplot(1,2,2);

surf(n_mesh,u_mesh,is_FS_overlap_K_dn_and_Kp_up(:,:)','edgecolor','none','FaceAlpha',0.5)
hold on
%contour(n_mesh,u_mesh,is_FS_overlap_K_dn_and_Kp_up(:,:)',[0.5,0.5],'linecolor','y','linewidth',2)
hold off
%shading interp
xlabel('$n_e [\rm{nm}^{-2}]$','interpreter','latex')
ylabel('$u [\rm{meV}]$','interpreter','latex')
view(2)
%mycolormap = customcolormap_preset('pink-white-green');
%colorbar;%('southoutside');
colormap(customcolormap_preset('red-yellow-blue'));%sp(2),
caxis_max = 1;
%clim([-caxis_max caxis_max])
grid off
ax = gca;
ax.FontSize = 20;

set(gcf, 'Position',  [100, 100, 600, 500])

% figure(3)
% sp(1) = subplot(1,2,1);
% surf(n_mesh,u_mesh,is_three_pocket_K_up(:,:)','edgecolor','none','FaceAlpha',0.6)
% hold off
% %shading interp
% xlabel('$n_e [\rm{nm}^{-2}]$','interpreter','latex')
% ylabel('$u [\rm{meV}]$','interpreter','latex')
% view(2)
% %mycolormap = customcolormap_preset('red-yellow-blue');
% colorbar;%('southoutside');
% colormap(sp(1),customcolormap_preset('red-yellow-blue'));
% caxis_max = 1;
% %clim([-caxis_max caxis_max])
% ax = gca;
% ax.FontSize = 20;
% 
% sp(2) = subplot(1,2,2);
% surf(n_mesh,u_mesh,is_three_pocket_K_up(:,:)','edgecolor','none','FaceAlpha',0.6)
% hold off
% %shading interp
% xlabel('$n_e [\rm{nm}^{-2}]$','interpreter','latex')
% ylabel('$u [\rm{meV}]$','interpreter','latex')
% view(2)
% colorbar;%('southoutside');
% colormap(sp(2),customcolormap_preset('orange-white-purple'));
% caxis_max = 1;
% %clim([-caxis_max caxis_max])
% ax = gca;
% ax.FontSize = 20;
% 
% set(gcf, 'Position',  [100, 100, 1000, 500])


figure(4)
% sp(1) = subplot(1,2,1);
% surf(n_mesh,u_mesh,eh_pockets_num_K_up(:,:)','edgecolor','none','FaceAlpha',0.5)
% hold off
% %shading interp
% xlabel('$n_e [\rm{nm}^{-2}]$','interpreter','latex')
% ylabel('$u [\rm{meV}]$','interpreter','latex')
% view(2)
% %mycolormap = customcolormap_preset('red-yellow-blue');
% colorbar;%('southoutside');
% colormap(sp(1),customcolormap_preset('red-yellow-blue'));
% %caxis_max = 1;
% %clim([-caxis_max caxis_max])
% grid off
% ax = gca;
% ax.FontSize = 20;
% 
% sp(2) = subplot(1,2,2);
surf(n_mesh,u_mesh,eh_pockets_num_K_dn(:,:)','edgecolor','none','FaceAlpha',0.5)
hold off
%shading interp
xlabel('$n_e [\rm{nm}^{-2}]$','interpreter','latex')
ylabel('$u [\rm{meV}]$','interpreter','latex')
view(2)
%colorbar;%('southoutside');
colormap(customcolormap_preset('red-yellow-blue'));%sp(2),
%caxis_max = 1;
%clim([-caxis_max caxis_max])
grid off
ax = gca;
ax.FontSize = 20;

set(gcf, 'Position',  [100, 100, 600, 500])


%%
u_measured = 65;%67 , 57 < <72
% 61 for lambda=3
% 57 for lambda=4
% 53 for lambda = 5 nkxy=100
% 49 for lambda = 6 nkxy=100
% 42.5 for lambda = 8 nkxy=60

D_measured = 1;
u_vs_D = u_measured/D_measured;

[num,txt,raw] = xlsread('extract_BBG_PIP2_Sym12.xls');
n_extract1=num(:,1);
D_extract1=num(:,2);
u_extract1 = u_vs_D*D_extract1;

[num,txt,raw] = xlsread('extract_BBG_SC.xls');
n_extract2=num(:,1);
D_extract2=num(:,2);
u_extract2 = u_vs_D*D_extract2;

[num,txt,raw] = xlsread('extract_BBG_Sym4_PIP2.xls');
n_extract3=num(:,1);
D_extract3=num(:,2);
u_extract3 = u_vs_D*D_extract3;

figure(2)
% subplot(1,2,1)
% hold on
% plot3(n_extract1,u_extract1,1000*ones(size(u_extract1)),...
%     'o','Color','k','MarkerSize',5,'MarkerFaceColor','none')
% plot3(n_extract2,u_extract2,1000*ones(size(u_extract2)),...
%     'o','Color','b','MarkerSize',5,'MarkerFaceColor','cyan')
% %yline(u_measured,':r')
% axis([min(n_list) max(n_list), min(u_list),max(u_list)])
% hold off
% subplot(1,2,2)
hold on
plot3(n_extract1,u_extract1,1000*ones(size(u_extract1)),...
    'o','Color','k','MarkerSize',8,'MarkerFaceColor','none')
plot3(n_extract2,u_extract2,1000*ones(size(u_extract2)),...
    'o','Color','b','MarkerSize',8,'MarkerFaceColor','cyan')
%plot3(n_extract3,u_extract3,1000*ones(size(u_extract2)),...
%    'o','Color','w','MarkerSize',8,'MarkerFaceColor','none')
%yline(u_measured,':r')
axis([min(n_list) max(n_list), min(u_list),max(u_list)])
hold off

figure(4)
% subplot(1,2,1)
% hold on
% plot3(n_extract1,u_extract1,1000*ones(size(u_extract1)),...
%     'o','Color','k','MarkerSize',5,'MarkerFaceColor','none')
% plot3(n_extract2,u_extract2,1000*ones(size(u_extract2)),...
%     'o','Color','b','MarkerSize',5,'MarkerFaceColor','cyan')
% %yline(u_measured,':r')
% axis([min(n_list) max(n_list), min(u_list),max(u_list)])
% hold off
% subplot(1,2,2)
hold on
plot3(n_extract1,u_extract1,1000*ones(size(u_extract1)),...
    'o','Color','k','MarkerSize',8,'MarkerFaceColor','none')
plot3(n_extract2,u_extract2,1000*ones(size(u_extract2)),...
    'o','Color','b','MarkerSize',8,'MarkerFaceColor','cyan')
plot3(n_extract3,u_extract3,1000*ones(size(u_extract2)),...
    'o','Color','w','MarkerSize',8,'MarkerFaceColor','none')
%scan_n = [-0.0105;-0.0075];
%scan_u = [u_measured;u_measured];
%plot3(scan_n,scan_u,1000*ones(size(scan_n)),'--k','linewidth',1) % scanned line
%plot3([-0.00997;-0.00965],u_measured*ones(2,1), 1000*ones(2,1),...
%    'o','Color','k','MarkerSize',8,'MarkerFaceColor',[0.7,0,0]) % two transition points on the line scanned

axis([min(n_list) max(n_list), min(u_list),max(u_list)])
hold off

%% plot polarization
% n_list = charge_sign*n_list;
% [n_mesh,D_mesh] = meshgrid(n_list,D_list);
% figure(2)
% surf(n_mesh,D_mesh,polarization(:,:)','edgecolor','none')
% hold off
% %shading interp
% xlabel('$n [\rm{nm}^{-2}]$','interpreter','latex')
% ylabel('$D [meV]$','interpreter','latex')
%
% view(2)
% mycolormap = customcolormap_preset('red-yellow-blue');
% colorbar;%('southoutside');
% colormap(mycolormap);
% caxis_max = 1;
% caxis([-caxis_max caxis_max])
%
% %zlabel('$\epsilon(k)$','interpreter','latex')

%--------------------------------------
% for i_band = target_band
% n(i_band) = sum(sum(epsilon(:,:,i_band) < mu))*dkx*dky/(4*pi^2);
% E_kinetic(i_band) = sum(sum(epsilon(:,:,i_band).*(epsilon(:,:,i_band) < mu)));
% end