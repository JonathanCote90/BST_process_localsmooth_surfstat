function varargout = process_localsmooth_surfstat( varargin )
% PROCESS_LOCALSMOOTH_SURFSTAT: Spatial smoothing of the sources using SurfStat (KJ Worsley).
% This process is vased on ssmotth_surfstat from Brainstorm by
% Authors: Peter Donhauser, Francois Tadel, 2015-2016
% Modified by Jonathan Cote 2020

% marco_method is a hidden variable with the following definition:
% 'if (nargin >= 1), [varargout{1:nargout}] = feval(varargin{1}, varargin{2:end});end'
% It basically runs the function with all the inputs if there`s more than
% one and won`t if there is no input.
eval(macro_method);
end

%% ===== GET DESCRIPTION =====
function sProcess = GetDescription()
% Description the process
sProcess.Comment     = 'Local smoothing';
sProcess.FileTag     = '';
sProcess.Category    = 'Filter';
sProcess.SubGroup    = 'Sources';
sProcess.Index       = 401;
sProcess.InputTypes  = {'results', 'timefreq'};
sProcess.OutputTypes = {'results', 'timefreq'};
sProcess.nInputs     = 1;
sProcess.nMinFiles   = 1;

% Definition of the options
% === SCOUTS
sProcess.options.scouts.Comment = '';
sProcess.options.scouts.Type    = 'scout';
sProcess.options.scouts.Value   = {};

% Default values for some options
sProcess.isSourceAbsolute = 1;
sProcess.processDim  = [];    % Do not split matrix


% Definition of the options
% === DESCRIPTION
sProcess.options.help.Comment = ['This process uses SurfStatSmooth (SurfStat, KJ Worsley).<BR><BR>' ...
    'The smoothing is based only on the surface topography, <BR>' ...
    'not the real geodesic distance between two vertices.<BR>', ...
    'The input in mm is converted to a number of edges based<BR>', ...
    'on the average distance between two vertices in the surface.<BR><BR>'];
sProcess.options.help.Type    = 'label';
% === FWHM (kernel size)
sProcess.options.fwhm.Comment = '<B>FWHM</B> (Full width at half maximum):  ';
sProcess.options.fwhm.Type    = 'value';
sProcess.options.fwhm.Value   = {3, 'mm', 0};
end

%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess)
Comment = 'Local smoothing';
end

%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs) %#ok<DEFNU>
    % Initialize returned variable
    OutputFiles = {};
    % Get scouts
    AtlasList = sProcess.options.scouts.Value;
    % Convert from older structure (keep for backward compatibility)
    if isstruct(AtlasList) && ~isempty(AtlasList)
        AtlasList = {'User scouts', {AtlasList.Label}};
    end
    % No scouts selected: exit
    if isempty(AtlasList) || ~iscell(AtlasList) || (size(AtlasList,2) < 2) || isempty(AtlasList{1,2})
        bst_report('Error', sProcess, [], 'No scout selected.');
        return;
    end
    
        % Get options
    FWHM = sProcess.options.fwhm.Value{1} / 1000;

    % ===== LOAD DATA =====
    % Load the surface filename from results file
    switch (file_gettype(sInputs.FileName))
        case {'results', 'link'}
            FileMat = in_bst_results(sInputs.FileName, 0, 'SurfaceFile', 'GridLoc', 'Atlas', 'nComponents', 'HeadModelType');
            nComponents = FileMat.nComponents;
        case 'timefreq'
            FileMat = in_bst_timefreq(sInputs.FileName, 0, 'SurfaceFile', 'GridLoc', 'Atlas', 'HeadModelType');
            nComponents = 1;
        otherwise
            error('Unsupported file format.');
    end
    % Error: cannot smooth results on volume grids
    if ~strcmpi(FileMat.HeadModelType, 'surface') % ~isempty(FileMat.GridLoc)
        bst_report('Error', sProcess, sInputs, 'Spatial smoothing is only supported for surface head models.');
        sInputs = [];
        return;
    % Error: cannot smooth results that are already based on atlases
    elseif ~isempty(FileMat.Atlas)
        bst_report('Error', sProcess, sInputs, 'Spatial smoothing is not supported for sources based on atlases.');
        sInputs = [];
        return;
    % Error: only for constrained sources
    elseif (nComponents ~= 1)
        bst_report('Error', sProcess, sInputs, ['This process is only available for source models with constrained orientations.' 10 'With unconstrained orientations, the source maps are usually already very smooth.']);
        sInputs = [];
        return;
    end
    
    
    [sSurf, iSurf] = bst_memory('LoadSurface', FileMat.SurfaceFile);
    actualIAtlas = find(strcmp({sSurf.Atlas.Name}, AtlasList{1,1}));
    sScoutsIndex = find(strcmp({sSurf.Atlas(actualIAtlas).Scouts.Label},AtlasList{1,2}));
    sScouts = sSurf.Atlas(actualIAtlas).Scouts(sScoutsIndex);
    
    %% Vertex removal copied segment from panel_scout
    
    
    iRemoveVert = setdiff(1:size(sSurf.Vertices,1),[sScouts.Vertices]);

    [Vertices, Faces, Atlas] = tess_remove_vert(sSurf.Vertices, sSurf.Faces, iRemoveVert, sSurf.Atlas);
    % Remove the handles of the scouts
    for iAtlas = 1:length(Atlas)
        for is = 1:length(Atlas(iAtlas).Scouts)
            Atlas(iAtlas).Scouts(is).Handles = [];
        end
%         if isfield(Atlas(iAtlas).Scouts, 'Handles');
%             Atlas(iAtlas).Scouts = rmfield(Atlas(iAtlas).Scouts, 'Handles');
%         end
    end
    % Build new surface
    tag = ['keep' sScouts.Label];
    tag = tag(find(~isspace(tag)));
    sSurfNew = db_template('surfacemat');
    sSurfNew.Comment  = [sSurf.Comment ' | ' tag];
    sSurfNew.Vertices = Vertices;
    sSurfNew.Faces    = Faces;
    sSurfNew.Atlas    = Atlas;
    sSurfNew.iAtlas   = sSurf.iAtlas;
    
        % === SAVE NEW FILE ===
    % Output filename
    NewTessFile = strrep(file_fullpath(sSurf.FileName), '.mat', ['_' tag '.mat']);
    NewTessFile = file_unique(NewTessFile);
    % Save file back
    bst_save(NewTessFile, sSurfNew, 'v7');
    % Get subject
    [sSubject, iSubject] = bst_get('SurfaceFile', sSurf.FileName);
    % Register this file in Brainstorm database
    db_add_surface(iSubject, NewTessFile, sSurfNew.Comment);
    % Re-open one to show the modifications
    %view_surface(NewTessFile);
    
    
    %% Project on new surface
    [projectedMap, ResultsMat] = bst_project_sources( file_fullpath(sInputs.FileName), NewTessFile, 0, 0 );
    
    %% Make Idx map
    outputDir = bst_fileparts(file_fullpath(sInputs.FileName));
    idxMap = ResultsMat;
    idxMap.ImageGridAmp = zeros(length(ResultsMat.ImageGridAmp),1);
    idxMap.Comment = [sInputs.Comment '_idxMap'];
    ProtocolInfo = bst_get('ProtocolInfo');
    %idxMap.SurfaceFile = file_win2unix([ProtocolInfo.SUBJECTS '\' sSurf.FileName]);
    idxMap.SurfaceFile = file_win2unix([sSurf.FileName]);
    idxMap.ImageGridAmp(1:length(sInputs.A(:,1)),1) = 1:length(sInputs.A(:,1));
    %idxMap.FileType = 'presults';
    %idxMap.A = (1:length(idxMap.A))';
    %db_add(iSubject,idxMap)
    currentFilename = sInputs.FileName(strfind(sInputs.FileName,'results'):end-16);
    participant = sInputs.FileName(1:strfind(sInputs.FileName,'results')-2);
    outputFilename = bst_process('GetNewFilename',[ProtocolInfo.STUDIES participant],currentFilename);
    bst_save(outputFilename,idxMap, 'v7')
    db_add_data(ProtocolInfo.iStudy, outputFilename, idxMap)
    
    %% Project Idx map
    [projectedIdxMap, ResultsIdxMat] = bst_project_sources(outputFilename,NewTessFile,0,0);
    
    %% Second half of smoothing, smoothing the projected map
    	% Load surface
    SurfaceMat = in_tess_bst(ResultsMat.SurfaceFile);    
      
    % ===== PROCESS =====
    % Convert surface to SurfStat format
    cortS.tri = SurfaceMat.Faces;
    cortS.coord = SurfaceMat.Vertices';

    % Get the average edge length
    [vi,vj] = find(SurfaceMat.VertConn);
    Vertices = SurfaceMat.VertConn;
    meanDist = mean(sqrt((Vertices(vi,1) - Vertices(vj,1)).^2 + (Vertices(vi,2) - Vertices(vj,2)).^2 + (Vertices(vi,3) - Vertices(vj,3)).^2));
    % FWHM in surfstat is in mesh units: Convert from millimeters to "edges"
    FWHMedge = FWHM ./ meanDist;
    
    % Display the result of this conversion
    msgInfo = ['Average distance between two vertices: ' num2str(round(meanDist*10000)/10) ' mm' 10 ...
               'SurfStatSmooth called with FWHM=' num2str(round(FWHMedge * 1000)/1000) ' edges'];
    bst_report('Info', sProcess, sInputs, msgInfo);
    disp(['SMOOTH> ' strrep(msgInfo, char(10), [10 'SMOOTH> '])]);   
    
    
    % Smooth surface
    for iFreq = 1:size(ResultsMat.ImageGridAmp,3)
        smoothedVals(:,:,iFreq) = SurfStatSmooth(ResultsMat.ImageGridAmp(:,:,iFreq)', cortS, FWHMedge)';
    end
    test =1;
    % Force the output comment
    %sInput.CommentTag = [sProcess.FileTag, num2str(FWHM*1000)];
    
    %% Returning on full size
    roundedChunkIdx = round(ResultsIdxMat.ImageGridAmp);
    idxChunkInt = unique(roundedChunkIdx);
    
    %testing field
    idxUniqueAndClosest = zeros(size(roundedChunkIdx));
    for rounded_counter = 1:length(roundedChunkIdx)
        howMany = sum(roundedChunkIdx == roundedChunkIdx(rounded_counter));
        if howMany == 1
            idxUniqueAndClosest(rounded_counter) = 1;
        else
           sameIdx = find(roundedChunkIdx == roundedChunkIdx(rounded_counter));
           [mini, mink] = min(abs(roundedChunkIdx(rounded_counter) - ResultsIdxMat.ImageGridAmp(sameIdx)));
           idxUniqueAndClosest(mink) = 1;
        end
    end
    
    %idxChunkIntBoul= ismember(ResultsIdxMat.ImageGridAmp,idxChunkInt);
    OutputFiles = sInputs;
    OutputFiles.A = zeros(size(sInputs.A));
    OutputFiles.A(roundedChunkIdx(logical(idxUniqueAndClosest))) = smoothedVals(logical(idxUniqueAndClosest));
    OutputFiles.CommentTag = [sProcess.FileTag, ['local smoothing ' num2str(FWHM*1000)]];
    
end 
    
    
%% Copied function
function [Vertices, Faces, Atlas] = tess_remove_vert(Vertices, Faces, iRemoveVert, Atlas)
% TESS_REMOVE_VERT: Remove some vertices from a tesselation
%
% Usage:  [Vertices, Faces, Atlas] = tess_remove_vert(Vertices, Faces, iRemoveVert, Atlas=[])
%
% INPUTS:
%     - Vertices    : [N,3] matrix
%     - Faces       : [M,3] matrix
%     - iRemoveVert : indices of vertices to remove
%     - Atlas       : Atlas structure
% @=============================================================================
% This function is part of the Brainstorm software:
% https://neuroimage.usc.edu/brainstorm
% 
% Copyright (c)2000-2020 University of Southern California & McGill University
% This software is distributed under the terms of the GNU General Public License
% as published by the Free Software Foundation. Further details on the GPLv3
% license can be found at http://www.gnu.org/copyleft/gpl.html.
% 
% FOR RESEARCH PURPOSES ONLY. THE SOFTWARE IS PROVIDED "AS IS," AND THE
% UNIVERSITY OF SOUTHERN CALIFORNIA AND ITS COLLABORATORS DO NOT MAKE ANY
% WARRANTY, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO WARRANTIES OF
% MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE, NOR DO THEY ASSUME ANY
% LIABILITY OR RESPONSIBILITY FOR THE USE OF THIS SOFTWARE.
%
% For more information type "brainstorm license" at command prompt.
% =============================================================================@
%
% Authors: Francois Tadel, 2008-2013
% Parse inputs
if (nargin < 4) || isempty(Atlas)
    Atlas = [];
end
% Re-numbering matrix
iKeptVert = setdiff(1:size(Vertices,1), iRemoveVert);
iVertMap = zeros(1, size(Vertices,1));
iVertMap(iKeptVert) = 1:length(iKeptVert);
% Remove vertices
Vertices(iRemoveVert,:) = [];
% Find the faces that contain removed vertices
iRemoveFace = find(sum(ismember(Faces, iRemoveVert), 2));
% Remove faces from list
Faces(iRemoveFace, :) = [];
% Renumber indices
Faces = iVertMap(Faces);
% === UPDATE ATLAS ===
% Loop on the atlases
for iAtlas = 1:length(Atlas)
    % Initialize list of scouts to remove for this atlas
    iScoutRm = [];
    % Loop on the scouts
    for iScout = 1:length(Atlas(iAtlas).Scouts)
        % Remove vertices
        Atlas(iAtlas).Scouts(iScout).Vertices = setdiff(Atlas(iAtlas).Scouts(iScout).Vertices, iRemoveVert);
        Atlas(iAtlas).Scouts(iScout).Seed     = setdiff(Atlas(iAtlas).Scouts(iScout).Seed, iRemoveVert);
        % Renumber remaining vertices
        Atlas(iAtlas).Scouts(iScout).Vertices = iVertMap(Atlas(iAtlas).Scouts(iScout).Vertices);
        % Remove scout if there are no vertices left
        if isempty(Atlas(iAtlas).Scouts(iScout).Vertices)
            iScoutRm = [iScoutRm, iScout];
        % Set a new seed if necessary
        elseif isempty(Atlas(iAtlas).Scouts(iScout).Seed)
            Atlas(iAtlas).Scouts(iScout) = panel_scout('SetScoutsSeed', Atlas(iAtlas).Scouts(iScout), Vertices);
        % Renumber the seed
        else
            Atlas(iAtlas).Scouts(iScout).Seed = iVertMap(Atlas(iAtlas).Scouts(iScout).Seed);
        end
    end
    % Remove the empty scouts
    if ~isempty(iScoutRm)
        Atlas(iAtlas).Scouts(iScoutRm) = [];
    end
end
end

function [OutputFiles, ResultsMat] = bst_project_sources( ResultsFile, destSurfFile, isAbsoluteValues, isInteractive )
% BST_PROJECT_SOURCES: Project source files on a different surface (currents or timefreq).
%
% USAGE:  OutputFiles = bst_project_sources( ResultsFile, DestSurfFile, isAbsoluteValues=0, isInteractive=1 )
% 
% INPUT:
%    - ResultsFile      : Relative path to sources file to reproject
%    - destSurfFile     : Relative path to destination surface file
%    - isAbsoluteValues : If 1, interpolate absolute values of the sources instead of relative values
%    - isInteractive    : If 1, displays questions and dialog boxes
%                         If 0, consider that it is running from the process interface

% @=============================================================================
% This function is part of the Brainstorm software:
% https://neuroimage.usc.edu/brainstorm
% 
% Copyright (c)2000-2020 University of Southern California & McGill University
% This software is distributed under the terms of the GNU General Public License
% as published by the Free Software Foundation. Further details on the GPLv3
% license can be found at http://www.gnu.org/copyleft/gpl.html.
% 
% FOR RESEARCH PURPOSES ONLY. THE SOFTWARE IS PROVIDED "AS IS," AND THE
% UNIVERSITY OF SOUTHERN CALIFORNIA AND ITS COLLABORATORS DO NOT MAKE ANY
% WARRANTY, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO WARRANTIES OF
% MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE, NOR DO THEY ASSUME ANY
% LIABILITY OR RESPONSIBILITY FOR THE USE OF THIS SOFTWARE.
%
% For more information type "brainstorm license" at command prompt.
% =============================================================================@
%
% Authors: Francois Tadel, 2010-2016


ResultsFile = {ResultsFile};
%destSurfFile = {destSurfFile};
%% ===== PARSE INPUTS ======
if (nargin < 4) || isempty(isInteractive)
    isInteractive = 1;
end
if (nargin < 3) || isempty(isAbsoluteValues)
    isAbsoluteValues = 0;
end

%% ===== GROUP BY SURFACES =====
% Group the results files to process by the surface on which they were computed
% Objective: Computing only once the transformation for all the files in the same group
ResultsGroups = {};
SurfaceGroups = {};
nGroup = 0;
OutputFiles = {};
errMsg = [];
HeadModelType = [];
% Display progress bar
% if isInteractive
%     bst_progress('start', 'Project sources', 'Initialization...');
% end
% Get file type: results or timefreq
isTimefreq = strcmpi(file_gettype(ResultsFile{1}), 'timefreq');
% For each sources file: get surface
for iRes = 1:length(ResultsFile)
    % Read surface file
    if isTimefreq
        ResMat = in_bst_timefreq(ResultsFile{iRes}, 0, 'SurfaceFile', 'DataType', 'DataFile', 'HeadModelType');
        % Check the data type: timefreq must be source/surface based, and no kernel-based file
        if ~strcmpi(ResMat.DataType, 'results')
            errMsg = 'Only cortical maps can be projected.';
            bst_report('Error', 'process_project_sources', ResultsFile{iRes}, errMsg);
            continue;
        elseif ~isempty(strfind(ResultsFile{iRes}, '_KERNEL_'))
            errMsg = 'Cannot re-project kernel-based time-frequency cortical maps.';
            bst_report('Error', 'process_project_sources', ResultsFile{iRes}, errMsg);
            continue;
        elseif isempty(ResMat.SurfaceFile) && ~isempty(ResMat.DataFile)
            ResAssocMat = in_bst_results(ResMat.DataFile, 0, 'SurfaceFile');
            ResMat.SurfaceFile = ResAssocMat.SurfaceFile;
        end
        % Load related source file
        if ~isempty(ResMat.DataFile) && isempty(ResMat.HeadModelType)
            [AssociateMat,AssociateFile] = in_bst_results(ResMat.DataFile, 0, 'HeadModelType');
            ResMat.HeadModelType = AssociateMat.HeadModelType;
        end
    else
        ResMat = in_bst_results(ResultsFile{iRes}, 0, 'SurfaceFile', 'HeadModelType');
    end
    
    % Default HeadModelType: surface
    if ~isfield(ResMat, 'HeadModelType') || isempty(ResMat.HeadModelType)
        ResMat.HeadModelType = 'surface';
    % Else : Check the type of grid (skip volume head models)
    elseif ismember(ResMat.HeadModelType, {'volume'})
        wrnMsg = ['To project source grids, see online tutorial "Group analysis: Subjects corergistration":' 10 ...
                  'https://neuroimage.usc.edu/brainstorm/Tutorials/CoregisterSubjects#Volume_source_models' 10 ...
                  'Skipping file: "' ResultsFile{iRes} '"...'];
        if isInteractive
            disp(wrnMsg);
        else
            bst_report('Error', 'process_project_sources', ResultsFile{iRes}, wrnMsg);
        end
        continue;
    end
    % Check head model type: must be the same for all the files
    if isempty(HeadModelType)
        HeadModelType = ResMat.HeadModelType;
    elseif ~isempty(ResMat.HeadModelType) && ~isequal(HeadModelType, ResMat.HeadModelType)
        errMsg = 'All the source files must be of the type (surface, volume or mixed).';
        bst_report('Error', 'process_project_sources', ResultsFile{iRes}, errMsg);
        continue;
    end
    
    % Associated surface not defined: error
    if isempty(ResMat.SurfaceFile)
        errMsg = 'Associated surface file is not defined.';
        bst_report('Error', 'process_project_sources', ResultsFile{iRes}, errMsg);
        continue;
    end
    % Check that it is not the destination surface
    if file_compare(ResMat.SurfaceFile, destSurfFile)
        if isInteractive
            disp(['BST> WARNING: Source and destination surfaces are the same for file: ' ResultsFile{iRes}]);
        else
            errMsg = 'Source and destination surfaces are the same.';
            bst_report('Error', 'process_project_sources', ResultsFile{iRes}, errMsg);
        end
        continue;
    end
    % Look for surface filename in SurfaceGroups
    iGroup = find(file_compare(ResMat.SurfaceFile, SurfaceGroups));
    % If does not exist yet: create it
    if isempty(iGroup)
        nGroup = nGroup + 1;
        SurfaceGroups{nGroup} = ResMat.SurfaceFile;
        ResultsGroups{nGroup} = ResultsFile(iRes);
    % Group exist: add file to the group
    else
        ResultsGroups{iGroup}{end+1} = ResultsFile{iRes};
    end
end
% If destination surface = source surface for all files
if (nGroup == 0)
    if isempty(errMsg)
        errMsg = ['Source and destination surfaces are the same for all the selected files.' 10 'Nothing to project...'];
    end   
    if isInteractive
        bst_error(errMsg, 'Project sources', 0);
        bst_progress('stop');
    else
        bst_report('Error', 'process_project_sources', ResultsFile, errMsg);
    end
    return;
end
% Get protocol folders
ProtocolInfo = bst_get('ProtocolInfo');


%% ===== PROJECT SOURCES =====
isStopWarped = [];
iUpdatedStudies = [];
% Process each surface group
for iGroup = 1:nGroup
    % ===== GET INTERPOLATION =====
    srcSurfFile = SurfaceGroups{iGroup};
    if isInteractive
        bst_progress('start', 'Project sources', 'Loading surfaces...');
    end
    % Compute interpolation
    [WmatSurf, sSrcSubj, sDestSubj, srcSurfMat, destSurfMat, isStopWarped] = tess_interp_tess2tess(srcSurfFile, destSurfFile, isInteractive, isStopWarped);
    % Source subject and destination subject are the same
    isSameSubject = file_compare(sSrcSubj.FileName, sDestSubj.FileName);

    % ===== CREATE GROUP ANALYSIS SUBJECT =====
    % If src and dest subjects are not the same: create a "group analysis" subject
    if ~isSameSubject
        % If the destination is the group analysis subject: replace its name
        if strcmpi(sSrcSubj.Name, bst_get('DirDefaultSubject'))
            sSrcSubj.Name = bst_get('NormalizedSubjectName');
        end
        % If the destination is the group analysis subject: replace its name
        if strcmpi(sDestSubj.Name, bst_get('DirDefaultSubject'))
            sDestSubj = bst_get('NormalizedSubject');
            sDestSubj.Name = bst_get('NormalizedSubjectName');
        end
    end
    
    % ===== PROCESS EACH FILE =====
    nFile = length(ResultsGroups{iGroup});
    if isInteractive
        bst_progress('start', 'Project sources', 'Projecting sources...', 0, nFile);
    end
    % Process each results file in group
    for iFile = 1:nFile
        % Progress bar
        ResultsFile = ResultsGroups{iGroup}{iFile};
        if isInteractive
            bst_progress('inc', 1);
        end
        bst_progress('text', sprintf('Processing file #%d/%d: %s', iFile, nFile, ResultsFile));
        
        % ===== OUTPUT STUDY =====
        % Get source study
        [sSrcStudy, iSrcStudy] = bst_get('AnyFile', ResultsFile);
        % If result has to be save in "group analysis" subject
        if ~isSameSubject
            % New condition name
            NewCondition = strrep(sSrcStudy.Condition{1}, '@raw', '');
            % Get condition
            [sDestStudy, iDestStudy] = bst_get('StudyWithCondition', [sDestSubj.Name '/' NewCondition]);
            % Create condition if doesnt exist
            if isempty(iDestStudy)
                iDestStudy = db_add_condition(sDestSubj.Name, NewCondition, 0);
                if isempty(iDestStudy)
                    error(['Cannot create condition: "' normSubjName '/' NewCondition '".']);
                end
                sDestStudy = bst_get('Study', iDestStudy);
            end
        % Else: use the source study as output study
        else
            sDestStudy = sSrcStudy;
            iDestStudy = iSrcStudy;
        end

        % ===== LOAD INPUT FILE =====
        % Time-freq files
        if isTimefreq
            % Load file
            ResultsMat = in_bst_timefreq(ResultsFile, 0);
            ResFile = ResultsFile;
            % Change number of sources
            ResultsMat.RowNames = 1:size(ResultsMat.TF,1);
        % Source files: Read full
        else
            % Load file
            [ResultsMat,ResFile] = in_bst_results(ResultsFile, 1);
            % If it depends on the data file: try to use the comment of the data file
            if ~isempty(ResultsMat.DataFile) && strcmpi(file_gettype(ResultsMat.DataFile), 'data')
                % Find parent file in database
                [sStudy, iStudy, iData] = bst_get('DataFile', ResultsMat.DataFile);
                % If file was found: use its comment
                if ~isempty(sStudy) && ~isempty(sStudy.Data(iData).Comment)
                    ResultsMat.Comment = sStudy.Data(iData).Comment;
                end
            end
            % Unconstrained sources: Make a flat map
            if (ResultsMat.nComponents ~= 1)
                % Display warning before flattening, as this is not obvious
                strWarn = 'Unconstrained source maps are flattened (norm of the three orientations) before projection.';
                bst_report('Warning', 'process_project_sources', ResultsFile, strWarn);
                disp(['BST> Warning: ' strWarn]);
                % Apply norm of the three orientations
                ResultsMat = process_source_flat('Compute', ResultsMat);
            % Compute absolute values
            elseif isAbsoluteValues
                ResultsMat.ImageGridAmp = abs(ResultsMat.ImageGridAmp);               
            end

            ResultsMat.ChannelFlag = [];
            ResultsMat.GoodChannel = [];
        end
        % Remove link with original file
        ResultsMat.DataFile = [];
        % Check if the file was reprojected on an atlas
        if isfield(ResultsMat, 'Atlas') && ~isempty(ResultsMat.Atlas)
            wrnMsg = ['Cannot process atlas-based source files: Skipping file "' ResultsFile '"...'];
            if isInteractive
                disp(wrnMsg);
            else
                bst_report('Error', 'process_project_sources', ResultsFile, wrnMsg);
            end
            continue;
        end

        % ===== INTERPOLATION MATRIX =====
        % Surface source model: Simply use the surface-surface interpolation
        if strcmpi(HeadModelType, 'surface')
            Wmat = WmatSurf;
        % Mixed source models: Compute volume grid interpolations
        elseif strcmpi(HeadModelType, 'mixed')
            % Load MRI files
            sMriSrc  = in_mri_bst(sSrcSubj.Anatomy(sSrcSubj.iAnatomy).FileName);
            sMriDest = in_mri_bst(sDestSubj.Anatomy(sDestSubj.iAnatomy).FileName);
            % Compute interpolation
            [Wmat, destGridAtlas, destGridLoc, destGridOrient] = tess_interp_mixed(ResultsMat, WmatSurf, srcSurfMat, destSurfMat, sMriSrc, sMriDest, isInteractive);
            % Update output structure
            ResultsMat.GridAtlas  = destGridAtlas;
            ResultsMat.GridLoc    = destGridLoc;
            ResultsMat.GridOrient = destGridOrient;
%             % Check if there is a "Source model" atlas available
%             iModelSrc  = find(strcmpi({srcSurfMat.Atlas.Name}, 'Source model'));
%             iModelDest = find(strcmpi({destSurfMat.Atlas.Name}, 'Source model'));
%             if isempty(iModelDest) && ~isempty(iModelSrc)
%                 destSurfMat.Atlas(end+1) = srcSurfMat.Atlas(iModelSrc);
%             end
        else
            error(['Unsupported head model type: ' HeadModelType]);
        end
        
        % ===== PROJECT SOURCE MAPS =====
        % Time-freq file
        if isTimefreq
            % Apply interpolation matrix
            tmpTF = zeros(size(Wmat,1), size(ResultsMat.TF,2), size(ResultsMat.TF,3));
            for iFreq = 1:size(ResultsMat.TF,3)
                tmpTF(:,:,iFreq) = Wmat * ResultsMat.TF(:,:,iFreq);
            end
            ResultsMat.TF = tmpTF;
            % PAC: Apply interpolation to all measures
            if isfield(ResultsMat, 'sPAC') && ~isempty(ResultsMat.sPAC)
                if isfield(ResultsMat.sPAC, 'NestingFreq') && ~isempty(ResultsMat.sPAC.NestingFreq)
                    ResultsMat.sPAC.NestingFreq = Wmat * ResultsMat.sPAC.NestingFreq;
                end
                if isfield(ResultsMat.sPAC, 'NestedFreq') && ~isempty(ResultsMat.sPAC.NestedFreq)
                    ResultsMat.sPAC.NestedFreq = Wmat * ResultsMat.sPAC.NestedFreq;
                end
                if isfield(ResultsMat.sPAC, 'DirectPAC') && ~isempty(ResultsMat.sPAC.DirectPAC)
                    tmpTF = zeros(size(Wmat,1), size(ResultsMat.sPAC.DirectPAC,2), size(ResultsMat.sPAC.DirectPAC,3), size(ResultsMat.sPAC.DirectPAC,4));
                    for iLow = 1:size(ResultsMat.sPAC.DirectPAC,3)
                        for iHigh = 1:size(ResultsMat.sPAC.DirectPAC,4)
                            tmpTF(:,:,iLow,iHigh) = Wmat * ResultsMat.sPAC.DirectPAC(:,:,iLow,iHigh);
                        end
                    end
                    ResultsMat.sPAC.DirectPAC = tmpTF;
                end
            end
            ResultsMat.HeadModelType = HeadModelType;
            if ~isempty(ResultsMat.GridLoc) && (length(ResultsMat.RowNames) ~= size(ResultsMat.GridLoc,1))
                ResultsMat.RowNames = 1:size(ResultsMat.GridLoc,1);
            elseif isequal(ResultsMat.DataType, 'results') && isnumeric(ResultsMat.RowNames) && (length(ResultsMat.RowNames) ~= size(ResultsMat.TF,1))
                ResultsMat.RowNames = 1:size(ResultsMat.TF,1);
            end
        % Source file
        else
            % Apply interpolation matrix
            ResultsMat.ImageGridAmp = muliplyInterp(Wmat, double(ResultsMat.ImageGridAmp), ResultsMat.nComponents);
            
            % Apply interpolation to standart deviation matrix
            if isfield(ResultsMat, 'Std') && ~isempty(ResultsMat.Std)
                ResultsMat.Std = muliplyInterp(Wmat, double(ResultsMat.Std), ResultsMat.nComponents);
            end
            ResultsMat.ImagingKernel = [];
        end
        
        % === SAVE NEW RESULTS ===
        % Get source filename
        [tmp__, oldBaseName] = bst_fileparts(ResFile);
        % Remove KERNEL tag for saving full source files
        if ~isempty(strfind(oldBaseName, '_KERNEL_')) && isfield(ResultsMat, 'ImageGridAmp') && ~isempty(ResultsMat.ImageGridAmp)
            oldBaseName = strrep(oldBaseName, '_KERNEL', '');
        end
        % Prepare structure to be saved
        ResultsMat.SurfaceFile = destSurfFile;
        if ~isSameSubject
            ResultsMat.Comment = [sSrcSubj.Name '/' ResultsMat.Comment];
            newResultsFile = sprintf('%s_%s.mat', oldBaseName, file_standardize(sSrcSubj.Name));
        else
            ResultsMat.Comment = [ResultsMat.Comment ' | ' destSurfMat.Comment];
            newResultsFile = sprintf('%s_%dV.mat', oldBaseName, length(destSurfMat.Vertices));
        end
        % Surface file
        ResultsMat.SurfaceFile = file_win2unix(destSurfFile);
        % History: project source
        ResultsMat = bst_history('add', ResultsMat, 'project', ['Project sources: ' srcSurfFile ' => ' ResultsMat.SurfaceFile]);
        % Build full filename
        newResultsFileFull = bst_fullfile(ProtocolInfo.STUDIES, bst_fileparts(sDestStudy.FileName), newResultsFile);    
        newResultsFileFull = file_unique(newResultsFileFull);
        newResultsFile = file_short(newResultsFileFull);
        % Save new results file
        bst_save(newResultsFileFull, ResultsMat, 'v6');

        % === ADD FILE IN DATABASE ===
        % Create Results/Timefreq structure for database
        if isTimefreq
            sNewResults = db_template('Timefreq');
            sNewResults.FileName = newResultsFile;
            sNewResults.Comment  = ResultsMat.Comment;
            sNewResults.DataFile = ResultsMat.DataFile;
            sNewResults.DataType = ResultsMat.DataType;
            % If filename already exists in this study
            iExistingRes = find(file_compare({sDestStudy.Timefreq.FileName}, newResultsFile));
            if ~isempty(iExistingRes)
                % Replace previous Results
                sDestStudy.Timefreq(iExistingRes) = sNewResults;
            else
                % Add new result
                sDestStudy.Timefreq(end + 1) = sNewResults;
            end
        else
            sNewResults = db_template('Results');
            sNewResults.FileName = newResultsFile;
            sNewResults.Comment  = ResultsMat.Comment;
            sNewResults.DataFile = ResultsMat.DataFile;
            sNewResults.isLink   = 0;
            sNewResults.HeadModelType = HeadModelType;
            % If filename already exists in this study
            iExistingRes = find(file_compare({sDestStudy.Result.FileName}, newResultsFile));
            if ~isempty(iExistingRes)
                % Replace previous Results
                sDestStudy.Result(iExistingRes) = sNewResults;
            else
                % Add new result
                sDestStudy.Result(end + 1) = sNewResults;
            end
        end
        % Update study in database
        bst_set('Study', iDestStudy, sDestStudy);
        iUpdatedStudies = [iUpdatedStudies, iDestStudy];
        % Add to list of returned files
        OutputFiles{end+1} = newResultsFile;
    end
end


%% ===== UDPATE DISPLAY =====
if isInteractive
    bst_progress('stop');
end
if isempty(OutputFiles)
    return;
end
% Update node
panel_protocols('UpdateNode', 'Study', unique(iUpdatedStudies));
% Select first output study
panel_protocols('SelectStudyNode', iUpdatedStudies(1));
% Select first output file
panel_protocols('SelectNode', [], OutputFiles{1});
% Save database
db_save();


end



%% ===== APPLY INTERPOLATION MATRIX =====
function B = muliplyInterp(W, A, nComponents)
    switch (nComponents)
        case {0, 1}
            B = double(W * A);
        case 2
            B = zeros(2 * size(W,1), size(A,2));
            B(1:2:end,:) = double(W * A(1:2:end,:));
            B(2:2:end,:) = double(W * A(2:2:end,:));
        case 3
            B = zeros(3 * size(W,1), size(A,2));
            B(1:3:end,:) = double(W * A(1:3:end,:));
            B(2:3:end,:) = double(W * A(2:3:end,:));
            B(3:3:end,:) = double(W * A(3:3:end,:));
    end
    B = double(B);
end