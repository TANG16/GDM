function map = run_gdm_experiment_nii(varargin)

if nargin==0
    printhelp()
    return
end

if( strcmp(varargin{1},'--help') || isempty(varargin))
    printhelp()
    return;
end

if( strcmp(varargin{1},'-h') || isempty(varargin) )
    printhelp()
    return
end

if( strcmp(varargin{1},'--version') || isempty(varargin) )
    fprintf('Version 1.0.\n')
    return
end

if( strcmp(varargin{1},'-v') || isempty(varargin) )
    fprintf('Version 1.0.\n')
    return
end

if( strcmp(varargin{1},'-u') || isempty(varargin) )
    fprintf(' EXAMPLE USE (in matlab) \n');
    fprintf(' run_gdm_experiment_nii(''-i'',''test.csv'',''-o'',''.'',''-a'',1,''-b'',1) \n');
    fprintf(' EXAMPLE USE (in command line) \n');
    fprintf(' run_gdm_experiment_nii -i test.csv -o . -a 1 -b 1 \n');
    return
end

if( strcmp(varargin{1},'--usage') || isempty(varargin) )
    fprintf(' EXAMPLE USE (in matlab) \n');
    fprintf(' run_gdm_experiment_nii(''-i'',''test.csv'',''-o'',''.'',''-a'',1,''-b'',1) \n');
    fprintf(' EXAMPLE USE (in command line) \n');
    fprintf(' run_gdm_experiment_nii -i test.csv -o . -a 1 -b 1 \n');
    return
end


% function returns GDM statistical maps and associated p-values computed
% by analytically estimating permutation testing
%
% INPUT
%
% REQUIRED
% [--input, -i] : .csv file containing full paths to input images. (REQUIRED)
%              We assume that the first column contains subject identifying
%              information; the second column contains the path to the
%              images, while the last column contains label information.
%              First line of the file should contain header information.

% [--outputDir, -o] : directory where the output from all folds will be saved (REQUIRED)
%
% OPTIONAL
%
%
% [--a, -a] : discriminative regularization parameter (positive scalar) (default a=0.1)
% [--b, -b] : generative regularization parameter (positive scalar) (default b=0.1)
% [--usage, -u]  Prints basic usage message.          
% [--help, -h]  Prints help information.
% [--version, -v]  Prints information about software version.

%
% OUTPUT:
% map = structure that stores GDM statistics and p-values for every group/covariate
% provided with the following structure
% map.stat = cell that stores GDM statistic for groups/covariate in same
% order (output as .nii.gz file as well)
% map.p = cell that stores GDM statistic for groups/covariate in same
% order (output as .nii.gz file as well)
%
%
% NOTE: to compile this function do
% mcc -m  run_gdm_experiment_nii -A [NIFTI toolbox directory]
% EXAMPLE USE (in matlab)
% run_gdm_experiment_nii('-i','test.csv','-o','.','-a',1,'-b',1)
% EXAMPLE USE (in command line)
% run_gdm_experiment_nii -i test.csv -o . -a 1 -b 1


if( sum(or(strcmpi(varargin,'--input'),strcmpi(varargin,'-i')))==1)
    niiCSV=varargin{find(or(strcmpi(varargin,'--input'),strcmp(varargin,'-i')))+1};
else
    error('run_gdm_experiment_csv:argChk','Please specify input csv file!');
end


if( sum(or(strcmpi(varargin,'--outputDir'),strcmpi(varargin,'-o')))==1)
    outputDir=varargin{find(or(strcmp(varargin,'--outputDir'),strcmp(varargin,'-o')))+1};
else
    error('run_gdm_experiment_csv:argChk','Please specify output directory!');
end

if( sum(or(strcmpi(varargin,'--a'),strcmpi(varargin,'-a')))==1)
    params.lam1=varargin{find(or(strcmpi(varargin,'--a'),strcmp(varargin,'-a')))+1};
else
    params.lam1=0.1;
end

if( sum(or(strcmpi(varargin,'--b'),strcmpi(varargin,'-b')))==1)
    params.lam2=varargin{find(or(strcmpi(varargin,'--b'),strcmp(varargin,'-b')))+1};
else
    params.lam2=0.1;
end

if( sum(or(strcmpi(varargin,'--verbose'),strcmpi(varargin,'-vo')))==1)
    params.vo=varargin{find(or(strcmpi(varargin,'--verbose'),strcmp(varargin,'-vo')))+1};
else
    params.vo=0;
end

% create output directory
if (~exist(outputDir,'dir'))
    [status,~,~] = mkdir(outputDir);
    if (status == 0)
        error('run_gdm_experiment_csv:argChk','Cannot create output directory!');
    end
end

params.lam1=input2num(params.lam1);
params.lam2=input2num(params.lam2);
params.vo=input2num(params.vo);


% confirm validity of optional input arguments
validateFcn_C = @(x) (x>0);
validateFcn_vo = @(x) (x==0) || (x == 1);

if(~validateFcn_C(params.lam1))
    error('run_gdm_experiment_csv:argChk','Discriminative regularization parameter should be positive!');
end

if(~validateFcn_C(params.lam2))
    error('run_gdm_experiment_csv:argChk','Generative regularization parameter should be positive!');
end

if(~validateFcn_vo(params.vo))
    error('run_midas_experiment_csv:argChk','VO parameter should be either 0 or 1!');
end

disp('Done');
disp('GDM runs with the following parameteres');
disp(['niiCSV: ' niiCSV]);
disp(['OutputDir: ' outputDir]);
disp(['a: ' num2str(params.lam1)]);
disp(['b: ' num2str(params.lam2)]);
disp(['vo: ' num2str(params.vo)]);

% csv with features
fname=niiCSV;
if (~exist(fname,'file'))
    error('run_gdm_experiment_csv:argChk','Input feature .csv file does not exist');
end


% input data
% assumption is that the first column contains IDs, the second paths to
% files and the last contains labels
disp('Loading feature images...');
input=readtable(fname);
ID=input{:,1};
fnames=input{:,2};
Y=input{:,3:end};
covariates=input.Properties.VariableNames(3:end);

% read the images
count = size(Y,1);
info = load_untouch_nii_gz(fnames{1});

img.dimx = info.hdr.dime.dim(2) ;
img.dimy = info.hdr.dime.dim(3) ;
img.dimz = info.hdr.dime.dim(4) ;

X = zeros(count,img.dimx*img.dimy*img.dimz);
template = load_untouch_nii_gz(fnames{1});
foreground=zeros(size(template));
for ii=1:count
    disp(['Loading image ' fnames{ii} ]);
    nii = load_untouch_nii_gz(fnames{ii});
    X(ii,:) = nii.img(:) ;
    foreground=foreground+abs(nii.img);
end
foreground=double(foreground>0);


% z-score imaging features
X=zscore(X);

for i=1:size(X,1)
    i
    data{i}=reshape(X(i,:),size(foreground));
end

% clear X

toremove=find(isnan(sum(Y,2)));
if ~isempty(toremove)
    disp('Removed the following subjects due to incomplete information:')
    disp(ID(toremove));
%     data(isnan(sum(Y,2)))=[];
    X(isnan(sum(Y,2)),:)=[];
    Y(isnan(sum(Y,2)),:)=[];
end
rng('shuffle')

[J,A0,W0,C]=FBmodel_dual(X(:,foreground==1),zscore(Y(:,1)),zscore(Y(:,2:end)),params);

map.J = nan(size(foreground));
map.stat = nan(size(foreground));
map.p = nan(size(foreground));
map.J(foreground==1) = J;
map.stat(foreground==1) = J./sqrt(sum(C.^2,2));
map.p = 2*normcdf(-abs(map.stat),0,1);

Jh = X(:,foreground==1)'*(X(:,foreground==1)*J);
Ch = X(:,foreground==1)'*(X(:,foreground==1)*C);
map.Jhaufe = nan(size(foreground));
map.stathaufe = nan(size(foreground));
map.phaufe = nan(size(foreground));
map.Jhaufe(foreground==1) = Jh;
map.stathaufe(foreground==1) = Jh./sqrt(sum(Ch.^2,2));
map.phaufe = 2*normcdf(-abs(map.stathaufe),0,1);


disp('Saving results...')
if(params.vo==0)
    save([outputDir '/GDM_results.mat'],'map');
else
    save([outputDir '/GDM_results.mat'],'map');
end
disp('Done')
end

function mat2nii(map,template,prefix)
%function to write matlab map to a nifti file
%Input:
% 1) map: matlab map
% 2)template: a file to borrow nifti header from, one that matches the map
% dimensions
% 3) prefix: prefix of output
%Output:
% prefix.nii
% Example usage:
%mat2nii(map{1}.stat{1},'/cbica/software/external/fsl/4.1.5/data/standard/MNI152lin_T1_2mm_brain.nii.gz','output')

tmp=load_untouch_nii_gz(template);
tmp.hdr.dime.datatype=16;
tmp.img=map;

save_untouch_nii(tmp,prefix);
gzip([prefix '.nii'])
delete([prefix '.nii'])
disp(['Saved as ' prefix '.nii.gz'])
end

function printhelp()

fprintf(' function returns GDM statistical maps and associated p-values computed\n')
fprintf(' by analytically estimating permutation testing\n')
fprintf('\n')
fprintf(' INPUT\n')
fprintf('\n')
fprintf(' REQUIRED\n')
fprintf(' [--input, -i] : .csv file containing full paths to input images. (REQUIRED)\n')
fprintf('              We assume that the first column contains subject identifying\n')
fprintf('              information; the second column contains the path to the\n')
fprintf('              images, while the last column contains label information.\n')
fprintf('              First line of the file should contain header information.\n')
fprintf('\n')
fprintf(' [--outputDir, -o] : directory where the output from all folds will be saved (REQUIRED)\n')
fprintf('\n')
fprintf(' OPTIONAL\n')
fprintf('\n')
fprintf('\n')
fprintf(' [--a, -a] : discriminative regularization parameter (positive scalar) (default a=0.1)\n')
fprintf(' [--b, -b] : generative regularization parameter (positive scalar) (default b=0.1)\n')
fprintf(' [--usage, -u]  Prints basic usage message.          \n')
fprintf(' [--help, -h]  Prints help information.\n')
fprintf(' [--version, -v]  Prints information about software version.\n')
fprintf('\n')
fprintf('\n')
fprintf(' OUTPUT:\n')
fprintf(' map = structure that stores GDM statistics and p-values for every group/covariate\n')
fprintf(' provided with the following structure\n')
fprintf(' map.stat = cell that stores GDM statistic for groups/covariate in same\n')
fprintf(' order (output as .nii.gz file as well)\n')
fprintf(' map.p = cell that stores GDM statistic for groups/covariate in same\n')
fprintf(' order (output as .nii.gz file as well)\n')
fprintf('\n')
fprintf('\n')
fprintf(' NOTE: to compile this function do\n')
fprintf(' mcc -m  run_gdm_experiment_nii -A [NIFTI toolbox directory]\n')
fprintf(' EXAMPLE USE (in matlab)\n')
fprintf(' run_gdm_experiment_nii(''-i'',''test.csv'',''-o'',''.'',''-a'',1,''-b'',1)\n')
fprintf(' EXAMPLE USE (in command line)\n')
fprintf(' run_gdm_experiment_nii -i test.csv -o . -a 1 -b 1\n')
end

function o=input2num(x)
if isnumeric(x)
    o=x;
else
    o = str2double(x);
end
end