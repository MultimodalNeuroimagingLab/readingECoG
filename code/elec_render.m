%% Start here to make figures
% Freesurfer subjects directory
FSsubjectsdir = '/Applications/freesurfer/7.2.0/subjects/';
 
% load mni305 pial
[Lmnipial_vert,Lmnipial_face] = read_surf(fullfile(FSsubjectsdir,'fsaverage','surf','lh.pial'));
[Rmnipial_vert,Rmnipial_face] = read_surf(fullfile(FSsubjectsdir,'fsaverage','surf','rh.pial'));
 
 
%% Plot figure with left pial with electrodes in mni space
 
%  
% figure
% % transform to gifti
% gl.faces = Rmnipial_face+1;
% gl.vertices = Rmnipial_vert;
% gl = gifti(gl);
%  
% % render
% tH = ieeg_RenderGifti(gl);
% 
% % rotate with light
% ieeg_viewLight(90, 0)
%  
%% YBA Sub-01 right

figure

% transform to gifti
gl.faces = Rmnipial_face+1;
gl.vertices = Rmnipial_vert;
gl = gifti(gl);

% render
tH = ieeg_RenderGifti(gl);  

% RPTO 1-8
%ieeg_viewLight(180, 180) RO

RPTO1 = [1.1, -73.5, 13.2];
RPTO2 = [2, -81, 13.5];
RPTO3 = [3.5, -89, 13.8];
RPTO4 = [4.7, -97, 14.5];
RPTO5 = [15.5, -102, 12.1];
RPTO6 = [26.9, -99, 7.2];
RPTO7 = [39.1, -90, 0.79];
RPTO8 = [46, -86, 1.4];

ieeg_elAdd([RPTO1; RPTO2; RPTO3; RPTO4; RPTO5; RPTO6; RPTO7; RPTO8], 'k', 20)
ieeg_elAdd([RPTO1; RPTO2; RPTO3; RPTO4; RPTO5; RPTO6; RPTO7; RPTO8], [.99 .99 .99], 15)

% RP 1-6, RF 1-6
% ieeg_viewLight(90, 0) RL

RF1 = [40.6, -3.8, 59.3];
RF2 = [55.9, -3.8, 46.7];
RF3 = [57.8, -6.4, 32.5];
RF4 = [67.8, -10.5, 22.7];
RF5 = [68.9, -15.6, 14.5];
RF6 = [69.7, -21, 4.6];

ieeg_elAdd([RF1; RF2; RF3; RF4; RF5; RF6],'k', 20)
ieeg_elAdd([RF1; RF2; RF3; RF4; RF5; RF6], [.99 .99 .99], 15)

RP1 = [39.3, -44.5, 57.6];
RP2 = [58, -37.8, 48];
RP3 = [63.5, -34.6, 38];
RP4 = [65, -32.6, 23];
RP5 = [69, -29.8, 10.7];
RP6 = [70.2, -29, 2.6];

ieeg_elAdd([RP1; RP2; RP3; RP4; RP5; RP6], 'k', 20)
ieeg_elAdd([RP1; RP2; RP3; RP4; RP5; RP6], [.99 .99 .99], 15)

% RAIT 1-4, RMIT 1-4, RPIT 1-4
% ieeg_viewLight(-90, -90) RI

RAIT1 = [24.5, 7.5, -44.5];
RAIT2 = [36.7, 3.6, -45];
RAIT3 = [46.7, -3.5, -45.5];
RAIT4 = [52.6, -11.3, -43]; 

ieeg_elAdd([RAIT1; RAIT2; RAIT3; RAIT4], 'k', 20)
ieeg_elAdd([RAIT1; RAIT2; RAIT3; RAIT4], [.99 .99 .99], 15)

RMIT1 = [20.1, -30.8, -20.6]; 
RMIT2 = [28.6, -31.2, -22.8]; 
RMIT3 = [40.6, -31.5, -24.9];
RMIT4 = [51, -30.7, -30.8];

ieeg_elAdd([RMIT1; RMIT2; RMIT3; RMIT4], 'k', 20)
ieeg_elAdd([RMIT1; RMIT2; RMIT3; RMIT4], [.99 .99 .99], 15)

RPIT1 = [13, -57.7, -6.8]; 
RPIT2 = [24, -56.7, -9.5]; 
RPIT3 = [34, -55.6, -18.8];
RPIT4 = [44, -54.3, -22.2];

ieeg_elAdd([RPIT1; RPIT2; RPIT3; RPIT4], 'k', 20)
ieeg_elAdd([RPIT1; RPIT2; RPIT3; RPIT4], [.99 .99 .99], 15)

% rotate with light - view 1
% ieeg_viewLight(-45, -30)
% figurewrite(sprintf('Subj%d_R_view%d',1,1),-1,[],'rendering');
% exportgraphics(gcf,[sprintf('rendering/Subj%d_R_view%d.eps',1,1)]);

% rotate with light - view 2
% ieeg_viewLight(90, 0)
% figurewrite(sprintf('Subj%d_R_view%d',1,2),-1,[],'rendering');
% exportgraphics(gcf,[sprintf('rendering/Subj%d_R_view%d.eps',1,2)]);

% rotate with light - view 3
% ieeg_viewLight(-45, 30)

% rotate with light - view 4
% ieeg_viewLight(30, -30)

%% Significant Electrodes

ieeg_elAdd(RMIT3, 'k', 40 + 10 * elec_weights{1}(6))
ieeg_elAdd(RMIT3, [.99 .01 .01], 30 + 10 * elec_weights{1}(6))

ieeg_elAdd(RPIT2, 'k', 40 + 10 * elec_weights{1}(7))
ieeg_elAdd(RPIT2, [.99 .01 .01], 30 + 10 * elec_weights{1}(7))

ieeg_elAdd(RPIT3, 'k', 40 + 10 * elec_weights{1}(8))
ieeg_elAdd(RPIT3, [.99 .01 .01], 30 + 10 * elec_weights{1}(8))

ieeg_elAdd(RPIT4, 'k', 40 + 10 * elec_weights{1}(9))
ieeg_elAdd(RPIT4, [.99 .01 .01], 30 + 10 * elec_weights{1}(9))

ieeg_elAdd(RPTO4, 'k', 40 + 10 * elec_weights{1}(10))
ieeg_elAdd(RPTO4, [.99 .01 .01], 30 + 10 * elec_weights{1}(10))

ieeg_elAdd(RPTO5, 'k', 40 + 10 * elec_weights{1}(11))
ieeg_elAdd(RPTO5, [.99 .01 .01], 30 + 10 * elec_weights{1}(11))

ieeg_elAdd(RPTO6, 'k', 40 + 10 * elec_weights{1}(12))
ieeg_elAdd(RPTO6, [.99 .01 .01], 30 + 10 * elec_weights{1}(12))

ieeg_elAdd(RPTO7, 'k', 40 + 10 * elec_weights{1}(13))
ieeg_elAdd(RPTO7, [.99 .01 .01], 30 + 10 * elec_weights{1}(13))

ieeg_elAdd(RPTO8, 'k', 40 + 10 * elec_weights{1}(14))
ieeg_elAdd(RPTO8, [.99 .01 .01], 30 + 10 * elec_weights{1}(14))

% ieeg_elAdd_sizable([RMIT3; RPIT2; RPIT3; RPIT4; RPTO4; RPTO5; RPTO6; RPTO7; RPTO8], 30 + 10 * elec_weights{1}(6:14), [.99 .01 .01])

% rotate with light - view 1
% ieeg_viewLight(-60, -30)
% figurewrite(sprintf('Subj%d_R_view%d_sign',1,1),-1,[],'rendering');
% exportgraphics(gcf,[sprintf('rendering/Subj%d_R_view%d_sign.eps',1,1)]);

% rotate with light - view 2
% ieeg_viewLight(90, 0)
% figurewrite(sprintf('Subj%d_R_view%d_sign',1,2),-1,[],'rendering');
% exportgraphics(gcf,[sprintf('rendering/Subj%d_R_view%d_sign.eps',1,2)]);

%% YBA Sub-01 left

figure

% transform to gifti
gl.faces = Lmnipial_face+1;
gl.vertices = Lmnipial_vert;
gl = gifti(gl);
 
% render
tH = ieeg_RenderGifti(gl);

% LPT 1-6
% ieeg_viewLight(-90, 0) LL

LPT1 = [-50.5, -78.6, 10.4];
LPT2 = [-57, -70, 3.5];
LPT3 = [-62.3, -61, -0.3];  
LPT4 = [-65.2, -53.6, -3];     
LPT5 = [-66, -45, -7.3];
LPT6 = [-67, -36.5, -9.7];

ieeg_elAdd([LPT1; LPT2; LPT3; LPT4; LPT5; LPT6], 'k', 40)
ieeg_elAdd([LPT1; LPT2; LPT3; LPT4; LPT5; LPT6], [.99 .99 .99], 30)

% LAT 1-4, LAIT 1-4, LPIT 1-4
% ieeg_viewLight(-90, -90) LI

LAT1 = [-31.7, 18.4, -41];
LAT2 = [-38.3, 18.6, -39];
LAT3 = [-43.8, 14.7, -37.5];
LAT4 = [-49.2, 7.9, -37.3]; 

ieeg_elAdd([LAT1; LAT2; LAT3; LAT4], 'k', 40)
ieeg_elAdd([LAT1; LAT2; LAT3; LAT4], [.99 .99 .99], 30)

LAIT1 = [-28.5, -6.3, -37.8]; 
LAIT2 = [-36.1, -3.6, -47.3]; 
LAIT3 = [-44.2, -5.7, -46];
LAIT4 = [-49.3, -9.8, -43.5];

ieeg_elAdd([LAIT1; LAIT2; LAIT3; LAIT4], 'k', 40)
ieeg_elAdd([LAIT1; LAIT2; LAIT3; LAIT4], [.99 .99 .99], 30)

LPIT1 = [-31.6, -86, -19.5]; 
LPIT2 = [-37.6, -78.3, -17]; 
LPIT3 = [-41.1, -70.4, -20.5];
LPIT4 = [-43.7, -60.3, -23.5];

ieeg_elAdd([LPIT1; LPIT2; LPIT3; LPIT4], 'k', 40)
ieeg_elAdd([LPIT1; LPIT2; LPIT3; LPIT4], [.99 .99 .99], 30)


% rotate with light - view 1
% ieeg_viewLight(45, -30);
% figurewrite(sprintf('Subj%d_L_view%d',1,1),-1,[],'rendering');
% exportgraphics(gcf,[sprintf('rendering/Subj%d_L_view%d.eps',1,1)]);

% rotate with light - view 2
% ieeg_viewLight(-90, 0)
% figurewrite(sprintf('Subj%d_L_view%d',1,2),-1,[],'rendering');
% exportgraphics(gcf,[sprintf('rendering/Subj%d_L_view%d.eps',1,2)]);

%% Significant Electrodes


ieeg_elAdd(LPIT1, 'k', 40 + 10 * elec_weights{1}(6))
ieeg_elAdd(LPIT1, [.99 .01 .01], 30 + 10 * elec_weights{1}(6))



ieeg_elAdd([LPIT1; LPIT2; LPIT3; LPT3; LPT4], [.99 .01 .01], 30+elec_weights{1}(6:14))

% rotate with light - view 1
% ieeg_viewLight(60, -30);
% figurewrite(sprintf('Subj%d_L_view%d_sign',1,1),-1,[],'rendering');
% exportgraphics(gcf,[sprintf('rendering/Subj%d_L_view%d_sign.eps',1,1)]);

% rotate with light - view 2
% ieeg_viewLight(-90, 0)
% figurewrite(sprintf('Subj%d_L_view%d_sign',1,2),-1,[],'rendering');
% exportgraphics(gcf,[sprintf('rendering/Subj%d_L_view%d_sign.eps',1,2)]);

%% YBD Sub-02 left

figure

% transform to gifti
gl.faces = Lmnipial_face+1;
gl.vertices = Lmnipial_vert;
gl = gifti(gl);

% render
tH = ieeg_RenderGifti(gl);

% LF 1-4, LP 1-6
% ieeg_viewLight(-90, 0) LL

LF1 = [-60.3, 0.6, 36.2];
LF2 = [-62.3, -2.5, 29.9];
LF3 = [-65.8, -8.1, 21.9];
LF4 = [-66.8, -12, 15.1];

ieeg_elAdd([LF1; LF2; LF3; LF4], 'k', 40)
ieeg_elAdd([LF1; LF2; LF3; LF4], [.99 .99 .99], 30)

LP1 = [-45.9, -42.6, 61.9];
LP2 = [-50.4, -43, 54.4];
LP3 = [-56.5, -44.1, 48];
LP4 = [-59.3, -46.2, 40.7];
LP5 = [-62.5, -48.4, 31.7];
LP6 = [-63.5, -50.2, 22.2];

ieeg_elAdd([LP1; LP2; LP3; LP4; LP5; LP6], 'k', 40)
ieeg_elAdd([LP1; LP2; LP3; LP4; LP5; LP6], [.99 .99 .99], 30)

% LAT 1-6, LAST 1-4, LD 1-8, LMST 1-4, LPST 1-6
% ieeg_viewLight(-180, -90) LI

LAT1 = [-24.6, 11.7, -43.5];
LAT2 = [-30.7, 17, -42];
LAT3 = [-37.5, 19.3, -38.5];
LAT4 = [-44.5, 15.6, -35.8];
LAT5 = [-47.5, 10.6, -34.7];
LAT6 = [-54, 5.1, -34.2];

ieeg_elAdd([LAT1; LAT2; LAT3; LAT4; LAT5; LAT6], 'k', 40)
ieeg_elAdd([LAT1; LAT2; LAT3; LAT4; LAT5; LAT6], [.99 .99 .99], 30)

LAST1 = [-28.6, 8.3, -45];
LAST2 = [-35.1, 6.6, -42];
LAST3 = [-40.7, 5.3, -39.5];
LAST4 = [-46.8, 4.3, -38];

ieeg_elAdd([LAST1; LAST2; LAST3; LAST4], 'k', 40)
ieeg_elAdd([LAST1; LAST2; LAST3; LAST4], [.99 .99 .99], 30)

LD1 = [-16.8, -14.3, -31];
LD2 = [-20.8, -12, -36.5];
LD3 = [-25.7, -9.8, -38];
LD4 = [-31, -7.6, -36.5];
LD5 = [-37.7, -4.6, -47.5];
LD6 = [-43.1, -2.3, -45.5];
LD7 = [-48.6, 0, -37];
LD8 = [-52.7, 1.5, -36];

ieeg_elAdd([LD1; LD2; LD3; LD4; LD5; LD6; LD7; LD8], 'k', 40)
ieeg_elAdd([LD1; LD2; LD3; LD4; LD5; LD6; LD7; LD8], [.99 .99 .99], 30)

LMST1 = [-38, -14.3, -35.5];
LMST2 = [-42, -11.7, -44.5];
LMST3 = [-47.9, -8.6, -45];
LMST4 = [-54.3, -5.3, -34];

ieeg_elAdd([LMST1; LMST2; LMST3; LMST4], 'k', 40)
ieeg_elAdd([LMST1; LMST2; LMST3; LMST4], [.99 .99 .99], 30)

LPST6 = [-52, -25, -35.5];
LPST5 = [-53 -33, -32.5];
LPST4 = [-53.5, -42, -30.5];
LPST3 = [-53.5, -48.5, -25.5];
LPST2 = [-53, -55, -22];
LPST1 = [-50.5, -61, -13.5];

ieeg_elAdd([LPST1; LPST2; LPST3; LPST4; LPST5; LPST6], 'k', 40)
ieeg_elAdd([LPST1; LPST2; LPST3; LPST4; LPST5; LPST6], [.99 .99 .99], 30)

% LTO 1-8
% ieeg_viewLight(90, 0) LL

LTO1 = [-37.4, -87, 2.4];
LTO2 = [-26.5, -92.4, 1.8];
LTO3 = [-22.8, -99, 4];
LTO4 = [-15.6, -103, 5.1];
LTO5 = [-5.1, -102.8, 6.7];
LTO6 = [-2.7, -95.2, 8.3];
LTO7 = [-0.1, -88, 10];
LTO8 = [-1.4, -81.8, 11.7];

ieeg_elAdd([LTO1; LTO2; LTO3; LTO4; LTO5; LTO6; LTO7; LTO8], 'k', 40)
ieeg_elAdd([LTO1; LTO2; LTO3; LTO4; LTO5; LTO6; LTO7; LTO8], [.99 .99 .99], 30)

% rotate with light - view 1
% ieeg_viewLight(45, -30);
% figurewrite(sprintf('Subj%d_L_view%d',2,1),-1,[],'rendering');
% exportgraphics(gcf,[sprintf('rendering/Subj%d_L_view%d.eps',2,1)]);

% rotate with light - view 2
% ieeg_viewLight(-90, 0)
% figurewrite(sprintf('Subj%d_L_view%d',2,2),-1,[],'rendering');
% exportgraphics(gcf,[sprintf('rendering/Subj%d_L_view%d.eps',2,2)]);

%% Significant Electrodes

ieeg_elAdd([LPST1; LPST2; LPST3; LPST4; LP4; LTO6], 'k', 40+elec_weights{2}(1:6))
ieeg_elAdd([LPST1; LPST2; LPST3; LPST4; LP4; LTO6], [.99 .01 .01], 30+elec_weights{2}(1:6))

% rotate with light - view 1
% ieeg_viewLight(60, -30);
% figurewrite(sprintf('Subj%d_L_view%d_sign',2,1),-1,[],'rendering');
% exportgraphics(gcf,[sprintf('rendering/Subj%d_L_view%d_sign.eps',2,1)]);

% rotate with light - view 2
% ieeg_viewLight(-90, 0)
% figurewrite(sprintf('Subj%d_L_view%d_sign',2,2),-1,[],'rendering');
% exportgraphics(gcf,[sprintf('rendering/Subj%d_L_view%d_sign.eps',2,2)]);

%%

figure

% transform to gifti
gl.faces = Rmnipial_face+1;
gl.vertices = Rmnipial_vert;
gl = gifti(gl);

% render
tH = ieeg_RenderGifti(gl);

% RF 1-4, RP 1-6
% ieeg_viewLight(90, 0) RL

RF1 = [43, -10.4, 62.9];
RF2 = [48, -6.6, 54.8];
RF3 = [55.5, -4.5, 46.7];
RF4 = [59, -2.4, 38.6];

ieeg_elAdd([RF1; RF2; RF3; RF4], 'k', 40)
ieeg_elAdd([RF1; RF2; RF3; RF4], [.99 .99 .99], 30)

RP1 = [39.5, -44.8, 58.7];
RP2 = [56.5, -42.3, 45.8];
RP3 = [64, -40.7, 35.6];
RP4 = [65.5, -39.1, 25.2];
RP5 = [68, -35.1, 15.6];
RP6 = [70, -30.7, 6.7];

ieeg_elAdd([RP1; RP2; RP3; RP4; RP5; RP6], 'k', 40)
ieeg_elAdd([RP1; RP2; RP3; RP4; RP5; RP6], [.99 .99 .99], 30)

% RAT 1-6, RAST 1-4, RD 1-8, RMST 1-4, RPST 1-6
% ieeg_viewLight(-180, -90) RI

RAT1 = [31, 18, -40];
RAT2 = [35, 20.5, -38.5];
RAT3 = [39, 21.5, -36];
RAT4 = [45.5, 19.5, -31];

ieeg_elAdd([RAT1; RAT2; RAT3; RAT4], 'k', 40)
ieeg_elAdd([RAT1; RAT2; RAT3; RAT4], [.99 .99 .99], 30)

RAST1 = [28.6, 8.3, -45];
RAST2 = [35.2, 8.2, -40];
RAST3 = [41, 6.8, -41.5];
RAST4 = [47.2, 5.5, -40.5];
RAST5 = [52, 4.1, -38];
RAST6 = [57, 1.8, -35];

ieeg_elAdd([RAST1; RAST2; RAST3; RAST4; RAST5; RAST6], 'k', 40)
ieeg_elAdd([RAST1; RAST2; RAST3; RAST4; RAST5; RAST6], [.99 .99 .99], 30)

RD1 = [20.8, -12, -36.5];
RD2 = [25.7, -9.8, -38];
RD3 = [31, -7.6, -39.5];
RD4 = [37.7, -4.6, -48.5];
RD5 = [34.5, -6, -46];
RD6 = [42.5, -2.5, -46.5];
RD7 = [46.5, -2, -43.5];
RD8 = [50.7, 0.5, -36];

ieeg_elAdd([RD1; RD2; RD3; RD4; RD5; RD6; RD7; RD8], 'k', 40)
ieeg_elAdd([RD1; RD2; RD3; RD4; RD5; RD6; RD7; RD8], [.99 .99 .99], 30)

RMST1 = [38, -14.3, -36];
RMST2 = [42, -11.7, -44.5];
RMST3 = [47.9, -8.6, -45];
RMST4 = [54.3, -5.3, -34];

ieeg_elAdd([RMST1; RMST2; RMST3; RMST4], 'k', 40)
ieeg_elAdd([RMST1; RMST2; RMST3; RMST4], [.99 .99 .99], 30)

RPST6 = [52, -25, -34.5];
RPST5 = [53 -33, -31];
RPST5 = [53.5, -42, -28];
RPST3 = [53.5, -48.5, -24.5];
RPST2 = [53, -55, -22];
RPST1 = [50.5, -61, -13.5];

ieeg_elAdd([RPST1; RPST2; RPST3; RPST4; RPST5; RPST6], 'k', 40)
ieeg_elAdd([RPST1; RPST2; RPST3; RPST4; RPST5; RPST6], [.99 .99 .99], 30)

% RTO 1-8
% ieeg_viewLight(90, 0) RL

RTO1 = [8.5, -83, 4.2];
RTO2 = [4.5, -90.5, 5.3];
RTO3 = [8.5, -98.3, 6];
RTO4 = [11, -103.5, 7];
RTO5 = [21, -102, 6];
RTO6 = [28, -97.3, 3];
RTO7 = [32.8, -92.6, 2];
RTO8 = [42, -88, 1];

ieeg_elAdd([RTO1; RTO2; RTO3; RTO4; RTO5; RTO6; RTO7; RTO8], 'k', 40)
ieeg_elAdd([RTO1; RTO2; RTO3; RTO4; RTO5; RTO6; RTO7; RTO8], [.99 .99 .99], 30)

% rotate with light - view 1
% ieeg_viewLight(-45, -30)
% figurewrite(sprintf('Subj%d_R_view%d',2,1),-1,[],'rendering');
% exportgraphics(gcf,[sprintf('rendering/Subj%d_R_view%d.eps',2,1)]);

% rotate with light - view 2
% ieeg_viewLight(90, 0)
% figurewrite(sprintf('Subj%d_R_view%d',2,2),-1,[],'rendering');
% exportgraphics(gcf,[sprintf('rendering/Subj%d_R_view%d.eps',2,2)]);

%% Significant Electrodes

ieeg_elAdd([RPST3; RPST5; RTO6; RTO7; RTO8], 'k', 40+elec_weights{2}(7:11))
ieeg_elAdd([RPST3; RPST5; RTO6; RTO7; RTO8], [.99 .01 .01], 30+elec_weights{2}(7:11))

% rotate with light - view 1
% ieeg_viewLight(-60, -30)
% figurewrite(sprintf('Subj%d_R_view%d_sign',2,1),-1,[],'rendering');
% exportgraphics(gcf,[sprintf('rendering/Subj%d_R_view%d_sign.eps',2,1)]);

% rotate with light - view 2
% ieeg_viewLight(90, 0)
% figurewrite(sprintf('Subj%d_R_view%d_sign',2,2),-1,[],'rendering');
% exportgraphics(gcf,[sprintf('rendering/Subj%d_R_view%d_sign.eps',2,2)]);
