%% Start here to make figures
% Freesurfer subjects directory
FSsubjectsdir = '/Applications/freesurfer/7.2.0/subjects/';
 
% load mni305 pial
[Lmnipial_vert,Lmnipial_face] = read_surf(fullfile(FSsubjectsdir,'fsaverage','surf','lh.pial'));
[Rmnipial_vert,Rmnipial_face] = read_surf(fullfile(FSsubjectsdir,'fsaverage','surf','rh.pial'));
 
 
%% Plot figure with left pial with electrodes in mni space
 
 
figure
% transform to gifti
gl.faces = Lmnipial_face+1;
gl.vertices = Lmnipial_vert;
gl = gifti(gl);
 
% render
tH = ieeg_RenderGifti(gl);

% rotate with light
ieeg_viewLight(-180, -90)
 
%% YBA Sub-01

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

ieeg_elAdd([RPTO1; RPTO2; RPTO3; RPTO4; RPTO5; RPTO6; RPTO7; RPTO8] ,'y',30)


% RP 1-6, RF 1-6
% ieeg_viewLight(90, 0) RL

RF1 = [40.6, 6.5, 59.3];
RF2 = [55.9, 1.4, 45.7];
RF3 = [60.8, -1.8, 34.2];
RF4 = [66.8, -7.7, 26];
RF5 = [68.9, -13.6, 17.5];
RF6 = [69.7, -21, 4.6];

ieeg_elAdd([RF1; RF2; RF3; RF4; RF5; RF6] ,'y',30)

RP1 = [39.3, -44.5, 57.6];
RP2 = [58, -37.8, 48];
RP3 = [63.5, -34.6, 38];
RP4 = [65, -35.6, 20];
RP5 = [69, -29.8, 10.7];
RP6 = [70.2, -29, 2.6];

ieeg_elAdd([RP1; RP2; RP3; RP4; RP5; RP6] ,'y',30)

% RAIT 1-4, RMIT 1-4, RPIT 1-4
% ieeg_viewLight(-90, -90) RI

RAIT1 = [24.5, 7.5, -45];
RAIT2 = [34.3, 4.4, -45.2];
RAIT3 = [41.7, 0, -46];
RAIT4 = [48.2, -5.9, -44.5]; 

ieeg_elAdd([RAIT1; RAIT2; RAIT3; RAIT4] ,'y',30)

RMIT1 = [32.5, -44.9, -23]; 
RMIT2 = [39.1, -44.6, -24]; 
RMIT3 = [48.4, -44.4, -21.5];
RMIT4 = [56.3, -44.6, -26];

ieeg_elAdd([RMIT1; RMIT2; RMIT3; RMIT4] ,'y',30)

RPIT1 = [31.5, -60.5, -18]; 
RPIT2 = [40.1, -59, -22]; 
RPIT3 = [49.8, -57.5, -18.5];
RPIT4 = [56.2, -58, -19.5];

ieeg_elAdd([RPIT1; RPIT2; RPIT3; RPIT4] ,'b',30)

% LPT 1-6
% ieeg_viewLight(-90, 0) LL

LPT1 = [-50.5, -78.6, 10.4];
LPT2 = [-42, -68, 3.5];
LPT3 = [-62.3, -61, -0.3];  
LPT4 = [-65.2, -53.6, -3];     
LPT5 = [-66, -45, -7.3];
LPT6 = [-67, -36.5, -9.7];

ieeg_elAdd([LPT1; LPT2; LPT3; LPT4; LPT5; LPT6] ,'y',30)

% LAT 1-4, LAIT 1-4, LPIT 1-4
% ieeg_viewLight(-90, -90) LI

LAT1 = [-31.7, 18.4, -41];
LAT2 = [-38.3, 18.6, -39];
LAT3 = [-43.8, 14.7, -37.5];
LAT4 = [-49.2, 7.9, -37.3]; 

ieeg_elAdd([LAT1; LAT2; LAT3; LAT4] ,'y',30)

LAIT1 = [-28.5, -6.3, -37.8]; 
LAIT2 = [-36.1, -3.6, -47.3]; 
LAIT3 = [-44.2, -5.7, -46];
LAIT4 = [-49.3, -9.8, -43.5];

ieeg_elAdd([LAIT1; LAIT2; LAIT3; LAIT4] ,'y',30)

LPIT1 = [-31.6, -86, -19.5]; 
LPIT2 = [-37.6, -78.3, -17]; 
LPIT3 = [-41.1, -70.4, -20.5];
LPIT4 = [-43.7, -60.3, -23.5];

ieeg_elAdd([LPIT1; LPIT2; LPIT3; LPIT4] ,'y',30)

%% YBD Sub-02

% LF 1-4, LP 1-6
% ieeg_viewLight(-90, 0) LL

LF1 = [-59, 3.8, 36.1];
LF2 = [-61, -0.9, 29];
LF3 = [-65, -7.2, 22.9];
LF4 = [-66.8, -12, 15.1];

ieeg_elAdd([LF1; LF2; LF3; LF4] ,'y',30)

LP1 = [-46, -34.5, 53];
LP2 = [-58, -40.5, 47.6];
LP3 = [-59, -43.8, 41.7];
LP4 = [-62, -46.1, 34.1];
LP5 = [-63.5, -48.2, 26.7];
LP6 = [-64.5, -49.8, 19.1];

ieeg_elAdd([LP1; LP2; LP3; LP4; LP5; LP6] ,'y',30)

% LAT 1-6, LAST 1-4, LD 1-8, LMST 1-4, LPST 1-6
% ieeg_viewLight(-180, -90) LI

LAT1 = [-24.6, 11.7, -43.5];
LAT2 = [-30.7, 17, -42];
LAT3 = [-37.5, 19.3, -38.5];
LAT4 = [-44.5, 15.6, -35.8];
LAT5 = [-47.5, 10.6, -34.7];
LAT6 = [-54, 5.1, -34.2];

ieeg_elAdd([LAT1; LAT2; LAT3; LAT4; LAT5; LAT6] ,'y',30)

LAST1 = [-28.6, 8.3, -45];
LAST2 = [-35.1, 6.6, -42];
LAST3 = [-40.7, 5.3, -39.5];
LAST4 = [-46.8, 4.3, -38];

ieeg_elAdd([LAST1; LAST2; LAST3; LAST4] ,'y',30)

LD1 = [-16.8, -14.3, -31];
LD2 = [-20.8, -12, -36.5];
LD3 = [-25.7, -9.8, -38];
LD4 = [-31, -7.6, -36.5];
LD5 = [-37.7, -4.6, -47.5];
LD6 = [-43.1, -2.3, -45.5];
LD7 = [-48.6, 0, -37];
LD8 = [-52.7, 1.5, -36];

ieeg_elAdd([LD1; LD2; LD3; LD4; LD5; LD6; LD7; LD8] ,'y',30)

LMST1 = [-38, -14.3, -35.5];
LMST2 = [-42, -11.7, -44.5];
LMST3 = [-47.9, -8.6, -45];
LMST4 = [-54.3, -5.3, -34];

ieeg_elAdd([LMST1; LMST2; LMST3; LMST4] ,'y',30)

LPST1 = [-52, -25, -35.5];
LPST2 = [-53 -33, -32.5];
LPST3 = [-53.5, -42, -30.5];
LPST4 = [-53.5, -48.5, -25.5];
LPST5 = [-53, -55, -22];
LPST6 = [-50.5, -61, -13.5];

ieeg_elAdd([LPST1; LPST2; LPST3; LPST4; LPST5; LPST6] ,'y',30)

% LTO 1-8
% ieeg_viewLight(90, 0) LL

LTO1 = [0, -83, 18.5];
LTO2 = [-2.5, -90.5, 18.5];
LTO3 = [-9, -94, 18];
LTO4 = [-13.7, -99.5, 17.7];
LTO5 = [-21, -99.5, 16.5];
LTO6 = [-26.8, -97.3, 14.9];
LTO7 = [-32.8, -92.6, 12.4];
LTO8 = [-38, -86.1, 9.5];

ieeg_elAdd([LTO1; LTO2; LTO3; LTO4; LTO5; LTO6; LTO7; LTO8] ,'y',30)


% RF 1-4, RP 1-6
% ieeg_viewLight(90, 0) RL

RF1 = [43, -10.4, 62.9];
RF2 = [48, -6.6, 54.8];
RF3 = [55.5, -4.5, 46.7];
RF4 = [59, -2.4, 38.6];

ieeg_elAdd([RF1; RF2; RF3; RF4] ,'y',30)

RP1 = [39.5, -44.8, 58.7];
RP2 = [56.5, -42.3, 45.8];
RP3 = [64, -40.7, 35.6];
RP4 = [65.5, -39.1, 25.2];
RP5 = [68, -35.1, 15.6];
RP6 = [70, -30.7, 6.7];

ieeg_elAdd([RP1; RP2; RP3; RP4; RP5; RP6] ,'y',30)

% RAT 1-6, RAST 1-4, RD 1-8, RMST 1-4, RPST 1-6
% ieeg_viewLight(-180, -90) RI

RAT1 = [31, 18, -40];
RAT2 = [35, 20.5, -38.5];
RAT3 = [39, 21.5, -36];
RAT4 = [45.5, 19.5, -31];

ieeg_elAdd([RAT1; RAT2; RAT3; RAT4] ,'y',30)

RAST1 = [28.6, 8.3, -45];
RAST2 = [35.2, 8.2, -40];
RAST3 = [41, 6.8, -41.5];
RAST4 = [47.2, 5.5, -40.5];
RAST5 = [52, 4.1, -38];
RAST6 = [57, 1.8, -35];

ieeg_elAdd([RAST1; RAST2; RAST3; RAST4; RAST5; RAST6] ,'y',30)

RD1 = [20.8, -12, -36.5];
RD2 = [25.7, -9.8, -38];
RD3 = [31, -7.6, -39.5];
RD4 = [37.7, -4.6, -48.5];
RD5 = [34.5, -6, -46];
RD6 = [42.5, -2.5, -46.5];
RD7 = [46.5, -2, -43.5];
RD8 = [50.7, 0.5, -36];

ieeg_elAdd([RD1; RD2; RD3; RD4; RD5; RD6; RD7; RD8] ,'y',30)

RMST1 = [38, -14.3, -36];
RMST2 = [42, -11.7, -44.5];
RMST3 = [47.9, -8.6, -45];
RMST4 = [54.3, -5.3, -34];

ieeg_elAdd([RMST1; RMST2; RMST3; RMST4] ,'y',30)

RPST1 = [52, -25, -34.5];
RPST2 = [53 -33, -31];
RPST3 = [53.5, -42, -28];
RPST4 = [53.5, -48.5, -24.5];
RPST5 = [53, -55, -22];
RPST6 = [50.5, -61, -13.5];

ieeg_elAdd([RPST1; RPST2; RPST3; RPST4; RPST5; RPST6] ,'y',30)

% RTO 1-8
% ieeg_viewLight(90, 0) RL

RTO1 = [1, -83, 18.5];
RTO2 = [3.5, -90.5, 18.5];
RTO3 = [9, -94, 18];
RTO4 = [13.7, -99.5, 17.7];
RTO5 = [21, -99.5, 16.5];
RTO6 = [26.8, -97.3, 14.9];
RTO7 = [32.8, -92.6, 12.4];
RTO8 = [38, -86.1, 9.5];

ieeg_elAdd([RTO1; RTO2; RTO3; RTO4; RTO5; RTO6; RTO7; RTO8] ,'y',30)
