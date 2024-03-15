classdef fig_params
    properties

        %alpha for significance test
        alpha_noise = 0.1;
        
        %Line Plot Options
        p_line_width = 2.5;
        p_smooth_window = 30; %best if set to block size (30)
        p_line_style = '-';
        p_marker_size = 5;
        p_mark_split_width = 1.5;
        p_mark_split_col = [0.65 0.65 0.65];
        p_sig_line_col = [0 1 1;1 1 0];
        p_noise_col = [0 0 0];
        p_population_col_success = [0.98 0.41 0.17];
        p_population_col_fail = [0.25 0.25 0.25];
        p_population_alpha = 0.25;
        p_annealing_col = [0.25 0.25 0.25];
        p_class_col = [0.25 0.25 0.25; 0.98 0.41 0.17];
        p_marker='none';
        p_MarkerNpnts=10; % every how many points we are putting a marker
        
        %  error bar
        shade_alpha=0.5; % alpha of the shaded errorbar
        ErrBarColor='r'; % color of the error bar if it is not specified
        
        % significance level plots
        sig_line_width=[2 4 6]; % corresponding to [0.05 0.01 0.001]
        STD_method='sem';
        CircularData=0; % are we using circular data?
        UseCos=0; % are we usding cosd to report angles in degree?
        
        %Scatter Plot Options
        point_color = [0.75 0.75 0.75];
        point_size = 14;
        point_marker_repitition = 'o'; %'o'
        point_repition_size = 16;
        
        %patch plots
        patch_alpha = 0.25;
        patch_alpha_min = 0.25;
        patch_alpha_max = 0.25;
        
        %quiver plot
        q_window = 5; %number of samples to smooth for plotting shift
        q_head_length = 3;
        q_head_width = 3;
        q_scale = 0.9 %scale (0 to 1) of the length of the line between points
        
        %Target Specific Colors
        target_colors = [0.98 0.41 0.17;0.3 0.2 1];
        target_marker_scale = 15; %how much larger make the target marker than the others
        target_learn_colors = [0.98 0.41 0.17;0.3 0.2 1]; %[0.75 0.1 0;0 0.1 0.75];
        
        %bar plots
        bar_face_alpha = 0.75;
        bar_edge_alpha = 1;
        bar_width=0.5;
        bar_facecol=[0 0 0];
        bar_edgecol=[0 0 0];
        
        %histograms
        histplot_calperc=1; % show percentage instead of count for histogram plots

        %polar plots
        polar_marker='o'
        polar_ThetaColor='blue';
        polar_RColor='red';   
        
        % image plots
        imageplotfunc='pcolor';
       
        %Global figure options
        font_size = 10;%14
        font_name = 'Arial';
        font_weight = 'normal';
        TitleFontSize=8;
        units = 'centimeters';%'pixels';
        line_width =2;
        font_angle='normal';
        SigStar_fontsize=2;
        Xticks=[];
        XticksLabels=[];
        Yticks=[];
        YticksLabels=[];
        PaperQualityFigs=0; % are we generating high quality figures
        % imagesc plots 
        image_colormap='parula'; %magma(128)
        caxis_limits='auto'; % set limits of color axis
        OriginLine=1;   % plot two lines cross the origin [0 0]
        OriginLine_width=2; % width of origin line 
        daspect=[]; %aspect ratio of the axis
        enforce_daspect=0; %  are we enforcing a aspect ratio
        % figure saving
        SaveEachFrame=0;  % should we save each frame  format as well
        SaveEachFrameFrmt='epsc'
        SaveEachSubplots=1; % are we saving each of the figure's subplots?
        CreateFolder4Fig=1; % are creating folder fro each figure 
        FigDPI='-r500'; % DPI for saving the figure
        % subplot numbers
        ThisSubplot=''; % pass a 3 item vector for desired suybplot for example [3 2 1];
        Sp=[]; % subplot 
        AppendTitles=0; % aree we adding this title to the tile already in the plot 
        % axis definitions
        axis_type='square';
        LegendTxt=''; % legend of the current text 
        AxesGrid='off'; % is the axis grid on
        LegendLoc='best';
        LegendFontSiz=9;
        % handling the data 
        NormalizebyMax=0; % are we normalizing the mean data in PlotMeanStd by their maximum 
        NormalizebyMean=0; % are we normalizing the mean data in PlotMeanStd by their mean 
        SmoothingMethod % method to smooth the data
        WidthSmoothing % width of smoothing 
        WidthSmoothingDim2 % width of smoothing in second dimension 
        SmoothTimeAxis=0; % are smoothing the time axis as well
        IsthisAxisTime=1; % is this axis a time axis 
        include_n=1; % include n when ploting 
        NPnts_SubtractBaseLine=10; % how many data points are we subtracting the baseline?
        SubtractBaseLine=0; % are subtracting the baseline
        ForceShaded=0 % are we forcing this shaded type for plotmeanstd
        % title
        ThisTitle=''; % title of this plot
        ThisXLabel=''; % Xlabel of this plot 
        ThisYLabel=''; % Xlabel of this plot 
        %table
        TabelSheet=1; % current sheet of the table
        % stats
        performtrend_stattest=0; % are we performing the trend statistical test on the data to see if there are trends 
        trend_stattest_alpha=0.05; % aplha for trend stat test
    end %properties

    properties (Access=private)
        ManData=ManipulateData;
    end
    
    methods
        function  obj=fig_params(varargin)
            if nargin~=0 % initialize vars
                obj=obj.ParseParams(varargin) ; %%Process optional inputs
            end
        end
        function  obj=ParseParams(obj,InputArgs)
            %Process optional inputs
            if mod(length(InputArgs), 2) ~= 0, error('Must pass key/value pairs for options.'); end
            for i = 1:2:length(InputArgs)
                try
                    obj.(InputArgs{i}) = InputArgs{i+1};
                catch
                    error('Couldn''t set option ''%s''.', InputArgs{2*i-1});
                end
            end
        end
        %% figure functions
        function FormatAxes(obj,ax_handle,varargin)
            global AnalysisOpts
            obj=obj.ParseParams(varargin);
            obj.font_size=10;
            set(ax_handle,'fontsize',obj.font_size,...
                'fontname',obj.font_name,...
                'fontweight',obj.font_weight,...
                'linewidth',obj.line_width,...
                'box','off','units',obj.units);
            AspectRatio=ax_handle.DataAspectRatio;

          %  axis tight
            axis(obj.axis_type);
            if ~isempty(obj.daspect)
                daspect(obj.daspect)
            elseif obj.enforce_daspect % put this back if you are analysing behavior
                 daspect(AspectRatio)
            end
%             % Change the font size of subplot titles            
%             titleHandle = get(ax_handle, 'Title');  % Get the title handle
%             set(titleHandle, 'FontSize', obj.TitleFontSize);  % Set the title font size
%             drawnow
        end
        function FormatAllAxesFig(obj,h,varargin) % formats all of the subplots of the figure
            global AnalysisOpts
            obj=obj.ParseParams(varargin);
            Children=h.Children;
            figure(h);
            obj.font_size=10;
            if size(Children,1)==0
                obj.FormatAxes(gca);
            else
                AxesInd=arrayfun(@(x) strcmp(Children(x).Type,'axes'),1:size(Children,1));
                for spn=find(AxesInd)
                    subplot(Children(spn));
                    obj.FormatAxes(gca);
                end
            end
            % try to set all of the fot size
            set(findall(gcf, 'Type', 'axes'), 'fontsize',obj.font_size,...
                'fontname',obj.font_name,...
                'fontweight',obj.font_weight);
            % Change the font size of subplot titles
            subplotHandles = findall(gcf, 'type', 'axes');
            for i = 1:length(subplotHandles)
                titleHandle = get(subplotHandles(i), 'Title');  % Get the title handle
                set(titleHandle, 'FontSize', obj.TitleFontSize);  % Set the title font size
            end
        end
        function FormatPolarAxes(obj,ax_handle)
            set(ax_handle,'fontsize',obj.font_size,...
                'fontname',obj.font_name,...
                'fontweight',obj.font_weight,...
                'linewidth',obj.line_width,...
                'box','off');
            axis tight
        end
        function FigureSizing(obj,handles,ax_position,fig_position)
            for i = 1:numel(handles)
                set(0, 'currentfigure', handles(i));
                set(gca,'units',obj.unitss,'position',ax_position);
                if ~isempty(fig_position)
                    set(handles(i),'units',obj.unitss,...
                        'position',fig_position);
                end
            end
        end
        function MakeAxesEqual(obj,handles)
            for i = 1:numel(handles)
                set(0, 'currentfigure', handles(i));
                temp{i} = gca;
            end
            linkaxes([temp{:}]);
        end
        function SaveFigs(obj,handles,filetype,name,savedir,fancyflag)
            saveCurFigs(handles,filetype,name,savedir,fancyflag)
            close(handles)
        end
        function col = GenDefaultColorArray(obj,array_size)
            col = NaN(array_size,3);
            for i = 1:array_size
                col(i,:) = obj.default_color;
            end
        end
        function SaveTable(obj,table_handle,save_dir,varargin)
            obj=obj.ParseParams(varargin);
            filename = [save_dir,'.xls'];
            writetable(table_handle,filename,'Sheet',obj.TabelSheet,'WriteRowNames',false,'WriteVariableNames',true)
        end
        function SaveFigSeries(obj,SaveFileName,save_dir,Figs,varargin)
            global AnalysisOpts
            obj=obj.ParseParams(varargin);
            if isempty(SaveFileName) % then everything is in savedir
                [save_dir,SaveFileName]=fileparts(save_dir);
                save_dir=[save_dir AnalysisOpts.FS];
            end
            if obj.CreateFolder4Fig
                if obj.SaveEachFrame % if we are saving each frame then create a folder for this imageand save everything inside
                    if ~strcmp(save_dir(end),AnalysisOpts.FS);save_dir=[save_dir AnalysisOpts.FS];end
                    save_dir=[save_dir SaveFileName  AnalysisOpts.FS];
                    mkdir([save_dir]);
                end
            end

            if iscell(Figs)
                ExistFig=cellfun(@isempty,Figs);
                ObjFigs=cellfun(@isobject,Figs); % is this a figure object
                OKFig=and(~ExistFig,ObjFigs);
                Figs=Figs(OKFig);
            else
                Figs=arrayfun(@(x) x,Figs,'UniformOutput',0);
            end
            
            if isempty(Figs)
                return
            end
            FileAlreadyExists=exist([save_dir SaveFileName '.tif'],'file');
            fn=0; % frame number
            for i = 1:numel(Figs)
                if iscell(Figs{i});Figs{i}=Figs{i}{1};end
                figure(Figs{i})
                set(gcf,'renderer','painters')    
                set(Figs{i},'Units', obj.units);
                pos = get(Figs{i},'Position');
                set(Figs{i},'PaperPositionMode','Auto','PaperUnits',obj.units,'PaperSize',[pos(3), pos(4)])
                obj.FormatAllAxesFig(Figs{i}); % format all of the fonts of the figure                                             
                if i==1 & ~FileAlreadyExists
                    saveas(Figs{i},[save_dir SaveFileName],'tif')
                    if obj.SaveEachFrame
                        fn=1;
                        if obj.PaperQualityFigs
                            export_fig([save_dir SaveFileName '_f' num2str(fn)],'-depsc',Figs{i},'-nocrop')
                        else
                           % saveas(Figs{i},[save_dir SaveFileName '_f' num2str(fn)],obj.SaveEachFrameFrmt)  
                            print([save_dir SaveFileName '_f' num2str(fn)], '-dsvg', '-painters', Figs{i},obj.FigDPI);
                            obj.SaveFigureSubplots(Figs{i},[SaveFileName '_f' num2str(fn)],save_dir);
                            savefig(Figs{i},[save_dir SaveFileName '_f' num2str(fn)]);
                        end
                    end
                end
                if ~isempty(Figs{i}) && i>1 | FileAlreadyExists
                    frame = getframe(Figs{i});
                    im{i} = frame2im(frame);
                    imwrite(im{i},[save_dir SaveFileName '.tif'],'WriteMode','append')
                    if obj.SaveEachFrame
                        fn=fn+1;
                        if obj.PaperQualityFigs
                            export_fig([save_dir SaveFileName '_f' num2str(fn)],'-depsc',Figs{i},'-nocrop')
                        else
                           % saveas(Figs{i},[save_dir SaveFileName '_f' num2str(fn)],obj.SaveEachFrameFrmt)
                            print([save_dir SaveFileName '_f' num2str(fn)], '-dsvg', '-painters', Figs{i},obj.FigDPI);
                            obj.SaveFigureSubplots(Figs{i},[SaveFileName '_f' num2str(fn)],save_dir);
                            savefig(Figs{i},[save_dir SaveFileName '_f' num2str(fn)]);
                        end
                    end
                 end
            end
            
        end
        
        function SaveCurrentFigs(obj,SaveFileName,save_dir,varargin) % saves all of the root figrues
            global AnalysisOpts
            obj=obj.ParseParams(varargin);
            if isempty(SaveFileName) % then everything is in savedir
                [save_dir,SaveFileName]=fileparts(save_dir);
                save_dir=[save_dir AnalysisOpts.FS];
            end
            CurrFigs=handle(sort(double(findobj(0,'type','figure'))));
            fn=0;
            for i = 1:numel(CurrFigs)
                figure(CurrFigs(i))
                set(gcf,'renderer','painters')
                if i==1
                    saveas(CurrFigs(i),[save_dir SaveFileName],'tif')
                    if obj.SaveEachFrame
                        fn=1;
                        % saveas(Figs{i},[save_dir SaveFileName '_f' num2str(fn)],'epsc')
                        %saveas(CurrFigs(i),[save_dir SaveFileName '_f' num2str(fn)],'svg')
                        print([save_dir SaveFileName '_f' num2str(fn)], '-dsvg', '-painters', CurrFigs(i),obj.FigDPI);
                        obj.SaveFigureSubplots(Figs{i},[SaveFileName '_f' num2str(fn)],save_dir);
                    end
                end
                if ~isempty(CurrFigs(i)) && i>1
                    frame = getframe(CurrFigs(i));
                    im{i} = frame2im(frame);
                    imwrite(im{i},[save_dir SaveFileName '.tif'],'WriteMode','append')
                    if obj.SaveEachFrame
                        fn=fn+1;
                        %saveas(CurrFigs(i),[save_dir SaveFileName '_f' num2str(fn)],'epsc')
                        print([save_dir SaveFileName '_f' num2str(fn)], '-dsvg', '-painters', CurrFigs(i),obj.FigDPI);
                        obj.SaveFigureSubplots(Figs{i},[SaveFileName '_f' num2str(fn)],save_dir);
                    end
                end
            end
            
        end
        function saveCurFigs(handles,format,fn_base,save_dir,fancyflag)
            
            %Camden MacDowell 2018
            % @Synposis Saves figures in handles as type and to a specified dir. If dir is empty,
            % saves to the current directory. To combine into a single pdf, save figs as
            % pdfs and then use function 'CombinedPDFSInFolder'
            %
            %
            % @Inputs
            % @handles (required)
            % list of figure handes where handles = [figurehandl1 figurehandle2, ...]
            %
            % @format (required)
            % Save format (e.g. '-dpdf' (or '-pdf' if fancyflag ==1) or any other format specificed in print fnct
            %
            % @fn_base (required)
            % Base for the file name including the file extension.
            % This way all generated figures have the same
            % suffix (easier to combined). Prefix for each file is 1 through numel(handels).
            %
            % @Savedir (optional)
            % Directory to save the figures. if empty (def) use the current directory.
            %
            % @fancyflag (optional)
            % If you have complex images that are having issues saving correctly using typically saveas
            % or if you want to save high quality images for publication purposes, set fancyflag to 1
            % this will use the 'export_fig' function which is slow but better for custom rendering
            
            if nargin <4
                save_dir = pwd;
            end
            if nargin<5
                fancyflag=0;
            end
            
            
            startDir = pwd; %so you return user to start dir after saving figs
            cd(save_dir);
            
            for cur_f = 1:numel(handles)
                %Save off the figure
                fig = handles(cur_f);
                set(fig,'Units', obj.units);
                pos = get(fig,'Position');
                set(fig,'PaperPositionMode','Auto','PaperUnits',obj.units,'PaperSize',[pos(3), pos(4)])
                filename = sprintf('%d_%s',cur_f,fn_base);
                if fancyflag
                    export_fig(filename,format,fig,'-nocrop')
                else
                    %Save the figure, with same suffix and iterate numbers as the prefix
                    print(fig,sprintf('%d_%s',cur_f,fn_base),format,'-r0')
                end
            end
            
            cd(startDir) %return user to start dir
            
        end %function
        function CollectFiguresandSave(obj,InputFileNames,OutputfileName,fmt,Pages) % loads current figure files and saves them all together in one plot 
            %@InputFileNames(Cell Array) file name of input files with their path. 
            %@OutputfileName(Char) file name of output files with their path.
            %@OutputExt(Char) extention of the output file 
            %@fmt format of the file 
            %@ Pages which page in the Tiff file we are looking at
            if isempty(fmt);fmt='tif';end
            if isempty(Pages);Pages=1;end
            NImg=length(InputFileNames); % number of images 
            [save_dir,SaveFileName]=fileparts(OutputfileName);
            mkdir(save_dir)
            k=1;
            for i=1:NImg
                if exist(InputFileNames{i},'file')
                    for Pg=Pages
                        im{k}=imread(InputFileNames{i},Pg);
                        if i==1 & Pg==Pages(1)
                            imwrite(im{k},[OutputfileName '.tif'],'tif')
                        else
                            imwrite(im{k},[OutputfileName '.tif'],'WriteMode','append')
                        end
                    end
                    k=k+1;
                else
                    fprintf(1,'\nFile %s does not exist ...',InputFileNames{i});
                end
            end
                  
        end
        function CombinePDFsInFolder(suffix,output_fn,deleteFlag)
            
            %Camden MacDowell 2018
            %Combines all pdfs in a folder with the designated suffix into one pdf
            % deleteflag : 0 def, don't delete original pdfs, 1, delete original pdfs.
            %
            % output_fn - String of the complete output file name (including path). Leave blank to have it
            % be labeled 'CombinedPDF'.
            %
            % NOTE: Suffix is needs to contain both a '*' and '.pdf' '*suffix.pdf.'
            
            if nargin <2
                output_fn = [];
            end
            if nargin <3
                deleteFlag = 0;
            end
            
            if isempty(output_fn)
                output_fn = [pwd filesep 'CombinedPDF.pdf'];
            end
            
            %grab the files of interest
            [file_list, ~] = GrabAnyFileType(suffix, 0);
            
            %Now order by their number %works of the assumption that it's
            %directory/#_name
            if numel(file_list) > 9  %It'll already be the correct order unless more than 9 long (this is where the name gets messed up)
                for i = 1:numel(file_list)
                    [~, name] = fileparts(file_list{i});
                    temp = regexp(name,'_');
                    temp = temp(1)-1; %only get the position before the first _ instance
                    Order(i) = str2num(name(1:temp));
                end
                [~,indx] = sort(Order,'ascend');
                file_list = file_list(indx);
            end
            
            
            %Check if the combined file already exists
            if isfile(output_fn)
                warning('Output file already exists, overwriting file')
                delete(output_fn)
            end
            
            %Create combined PDF
            append_pdfs(output_fn, file_list{:});
            
            %delete individual files if desired.
            if deleteFlag
                for cur_f = 1:numel(file_list)
                    delete(file_list{cur_f});
                end
            end
            
            
        end %function end
        function SaveFigureSubplots(obj,h,SaveFileName,save_dir) 
            if ~obj.SaveEachSubplots;return;end
            try
                % make a new directory for the subplots if it doesn't exist
                save_dirSp=[save_dir 'Subplots'];
                if ~exist(save_dirSp,'file')
                    mkdir(save_dir,'Subplots')
                end
                
                
                % Find all subplot handles
                subplotHandles = findall(gcf, 'type', 'axes');
                
                % Set the magnification factor and common size (adjust as needed)
                magnificationFactor = 0.5;  % Increase the size by a factor of 2
                %  commonSize = [800, 600];  % Use your desired common size
                screenSize = get(0, 'ScreenSize');
                commonSize = [screenSize(3), screenSize(4)];
                % Save each subplot as a separate file
                for i = 1:length(subplotHandles)
                    subplotHandle = subplotHandles(i);
                    
                    % Set the figure units to pixels
                    fig = figure('Units', 'pixels', 'Position', [0, 0, commonSize * (magnificationFactor+0.3)]);
                    
                    % Copy the subplot to the new figure
                    newSubplot = copyobj(subplotHandle, fig);
                    
                    set(newSubplot,'Units','pixels')
                    % Adjust the position of the copied subplot within the new figure
                    set(newSubplot, 'Position', [200, 200, commonSize * magnificationFactor]); % Adjust these values as needed
                    
                    % adjust the font sizes again
                    set(findall(gcf, 'Type', 'axes'), 'fontsize',46,...
                        'fontname',obj.font_name,...
                        'fontweight',obj.font_weight);
                    
                    subplotHandlesTitle = findall(gcf, 'type', 'axes');
                    for ii = 1:length(subplotHandlesTitle)
                        titleHandle = get(subplotHandlesTitle(ii), 'Title');  % Get the title handle
                        set(titleHandle, 'FontSize', obj.TitleFontSize);  % Set the title font size
                    end
                    
                    % Set the filename for saving (you can adjust the format and name)
                    filename = [SaveFileName '_Sp', num2str(i)];
                    %
                    %                 % Save the current figure to the specified file
                    print([save_dirSp filesep filename], '-dsvg', '-painters', fig,obj.FigDPI);
                    %
                    % Close the temporary figure
                    close(fig);
                end
                
            catch me
                disp(me.message)
            end           
        end
        function ShowMarkers(obj,Markers,Text)
            
            arrayfun(@(x) plot(0,x,'Marker',Markers{x},'MarkerSize',20),1:length(Text));
            arrayfun(@(x) text(0.5,x,Text{x},'FontSize',20),1:length(Text));
            axis off
        end
        function CloseFigs(obj,Figs)
            ExistFig=cellfun(@isempty,Figs);
            Figs=Figs(~ExistFig);
            for i = 1:numel(Figs)
                if iscell(Figs{i});Figs{i}=Figs{i}{1};end
            end
            cellfun(@(x) close(x),Figs,'UniformOutput',0);
        end
        function CloseAllFigs(~) % closes all of the figures 
            delete(get(0,'Children'));
        end
        function FigHndl=RenderFigure(obj,NFig,Position)  % generates a fig and returns handle
            if NFig==0;FigHndl=[];return;end
            if isempty(Position)
                Position=[0 0 0.9 0.9];
            end
            for i=1:NFig
                FigHndl{i}=figure;
                set(FigHndl{i},'Units','Normalized','Position',Position)
            end
        end
        function [h,Sp]=RenderSubplots(obj,NRow,NCol,h,N)
            % N is the total number of subplots. 
            if isempty(h);h=obj.RenderFigure(1,[]);else;figure(h);end
            % if we have given one NCol or NRow and N
            if ~isempty(NRow) && isempty(NCol) 
                NCol=ceil(N/NRow);
            elseif isempty(NRow) && ~isempty(NCol)
                NRow=ceil(N/NCol);
            end           
            if isempty(NRow);NRow=0;end;if isempty(NCol);NCol=0;end                        
            if ~isempty(N) && N>NRow*NCol
                NRow=ceil(N/3);NCol=ceil(N/NRow);
            end          
            Sp=arrayfun(@(x) subplot(NRow,NCol,x),1:NRow*NCol);
        end
        function [kColors,ThisColor] = getColorPalet(~,K,varargin) %Grabs K colors from predetermined color palet
            global AnalysisOpts

            color_palet = AnalysisOpts.CurrColorPalett;

            color_palet = repmat(color_palet, ceil(K/size(color_palet,1)),1);
            kColors = color_palet(1:K,:);
            % is there are two inputs, then get a single color from the
            % palet
            if ~isempty(varargin);ThisColor=color_palet(varargin{1},:);end
        end
        function Col=getOthercolormap(~,N,M,FlipFlag) % get a color map form othercolor file https://www.mathworks.com/matlabcentral/fileexchange/30564-othercolor
            % N name of the color map in othercolor.m
            % M number of colors 
            % FlipFlag flips the order of the colors in the colormap 
            Col=othercolor(N,M);
            if FlipFlag;Col=Col(end:-1:1,:);end
        end
        function outCol=getSingleColor(obj,Col) % gets single color 4 this condition
            if isempty(Col) 
                outCol=obj.getColorPalet(1);
            elseif (isinteger(Col) | isfloat(Col)) & length(Col)==1 
                outCol=obj.getColorPalet(Col);
                outCol=outCol(Col,:);
            elseif strcmp(Col,'random')
                ColInd=randi(15);
                outCol=obj.getColorPalet(ColInd);
                outCol=outCol(ColInd,:);
            elseif strcmp(Col,'none')
                outCol='none';
            else
                outCol=Col;
            end
        end
        function [outLineStyle]=getSingleLineStyle(obj,LineStyle) % get single line style          
             if isempty(LineStyle) 
                outLineStyle=obj.getLineStylePalet(1);
            elseif (isinteger(LineStyle) | isfloat(LineStyle)) & length(LineStyle)==1 
                outLineStyle=obj.getLineStylePalet(LineStyle);
                outLineStyle=outLineStyle(LineStyle);
            elseif strcmp(LineStyle,'random')
                LineStyleInd=randi(4);
                outLineStyle=obj.getLineStylePalet(LineStyleInd);
                outLineStyle=outLineStyle(LineStyleInd);
            elseif strcmp(LineStyle,'none')
                outLineStyle='none';
            else
                outLineStyle=LineStyle;
             end
             if iscell(outLineStyle)
                 outLineStyle=outLineStyle{1};
             end
        end
        function [kLineStyle,ThisLineStyle] = getLineStylePalet(~,K,varargin) %Grabs K linestyles from predetermined linestyle palet
            
            linestyle_palet ={ '-', '--',':','-.'};
            linestyle_palet = repmat(linestyle_palet,1,ceil(K/size(linestyle_palet,2)));
            kLineStyle = linestyle_palet(1:K);
            % is there are two inputs, then get a single linestyle from the
            % palet
            if ~isempty(varargin);ThisLineStyle=linestyle_palet(varargin{1});end
        end
        function [kMarkers,ThisMarker] = getMarkerPalet(~,K,varargin) %Grabs K markers from predetermined marker palet
            
            marker_palet ={'o','+','*','x','s','d','^','v','p','>','h','<','.'};
            marker_palet = repmat(marker_palet,1,ceil(K/size(marker_palet,2)));
            kMarkers = marker_palet(1:K);
            % is there are two inputs, then get a single marker from the
            % palet
            if ~isempty(varargin);ThisMarker=marker_palet(varargin{1});end
        end
       
        function [Title,xlbl,ylbl]=DetermineTitle(obj,Title,xlbl,ylbl)
            if ~isempty(obj.ThisTitle)
                Title=obj.ThisTitle;
            end
            Title=obj.ManData.RepDashSpace(Title);
            if ~isempty(obj.ThisXLabel)
                xlbl=obj.ThisXLabel;
            end
            xlbl=obj.ManData.RepDashSpace(xlbl);
            if ~isempty(obj.ThisYLabel)
                ylbl=obj.ThisYLabel;
            end
            ylbl=obj.ManData.RepDashSpace(ylbl); 
            if isempty(xlbl);xlbl='';end
            if isempty(ylbl);ylbl='';end              
        end
        %% Plotting functions
        function [hp,Y,YSTD,Yorg,Ysmoothed,Xsmoothed,Title]=PlotMeanStd(obj,X,Y,YSTD,Xlbl,Ylbl,Col,Shaded,Title,varargin)
            global AnalysisOpts
            % plots mean and std of data. if YSTD is not passed it
            % calculates it 
            if isempty(obj.p_marker);obj.p_marker='none';end %inistialze marker
            if ~isempty(varargin) & numel(varargin)==1
                obj.p_line_style=varargin{1};
            else
                obj.p_line_style='-';                
                obj=obj.ParseParams(varargin) ; %%Process optional inputs
            end
            
            [Title,Xlbl,Ylbl]=obj.DetermineTitle(Title,Xlbl,Ylbl);
            
            n=size(Y,1); % number of samples
           
            %if isnumeric(Col) && length(Col)==1;[~,Col]=obj.getColorPalet(100,Col);end
            Col=obj.getSingleColor(Col);
            if ~isempty(obj.ThisSubplot);subplot(obj.ThisSubplot(1),obj.ThisSubplot(2),obj.ThisSubplot(3));end
            if ~isempty(obj.Sp);subplot(obj.Sp);end
            
            % if we are smoothing the data            
            if ~isempty(obj.SmoothingMethod) & obj.CircularData==0
                Y=obj.ManData.SmoothData(Y,obj.WidthSmoothing,'SmoothingMethod',obj.SmoothingMethod);
                if obj.SmoothTimeAxis & obj.IsthisAxisTime
                    % smooth X axis as well
                    X=obj.ManData.SmoothData(X,obj.WidthSmoothing,'SmoothingMethod',obj.SmoothingMethod);
                    warning('X axis might not be time but it has been smoothed');
                end
            elseif obj.CircularData==1
                warning('Circular data has not been smoothed');
            end
            Ysmoothed=Y;Xsmoothed=X; % take a snapshot of the smoothed data as well
            
            % normalize the data 
            [Y,obj]=obj.ManData.NormalizeData(X,Y,obj);
            % calculate mean and STD               
            [Yorg,Y,YSTD]=obj.ManData.CalDataSTD(Y,YSTD,obj);
            
            if isempty(X)
                Y=Y(~isnan(Y));
                YSTD=YSTD(~isnan(Y));
                X=1:length(Y);
            end
            % change time axis if we need to 
            TimeInd=obj.LimitTimeAxis(X,Xlbl); X=X(TimeInd);Y=Y(TimeInd);YSTD=YSTD(TimeInd);
           
            % add calculations of latency 
            if obj.NormalizebyMax | obj.NormalizebyMean
                [OnSetLatencySec,PeakLatencySec]=obj.ManData.CalPeakBaselineLatency(Y,X); 
                obj.LegendTxt=[obj.LegendTxt sprintf(' %0.3f-%0.3f',OnSetLatencySec,PeakLatencySec)];
            end
            
            % if we are performing trend statistical test run it and add results to the title 
            if obj.performtrend_stattest
                TitleTrendStatTest=obj.ManData.PerformAllTrend_StatTest(X,Y,obj.trend_stattest_alpha);
                Title=[Title;TitleTrendStatTest];
            end
            
            % before adding get the active legends in case we want to add
            % legends in the future 
            [ActObj,ActObjLegends]=obj.GetActiveLegends([]);
            if length(X)==1 & ~obj.ForceShaded;Shaded=0;end
            if Shaded==1 % shaded
                shadedErrorBar(X,Y,YSTD,'lineProps',{'color',Col,'LineStyle',obj.p_line_style,...
                    'linewidth',obj.p_line_width},'patchSaturation',0.2); 
                alpha(obj.patch_alpha);
                hold on;
                hp=plot(X, Y, 'LineStyle',obj.p_line_style, 'Color', Col,'linewidth',obj.p_line_width); 
                hp.Visible='off';
                 % now add a merker every N Trials
                data_counter = 0;
                % Plot the data points with shading and markers at specific intervals
                for i = 1:length(X)-1
                    data_counter = data_counter + 1;
                    if data_counter == obj.p_MarkerNpnts
                        % Add a marker
                        plot(X(i), Y(i), 'Color', Col,'Marker',obj.p_marker,'MarkerSize',obj.p_marker_size);
                        data_counter = 0;  % Reset the counter
                    end
                end
            elseif Shaded==0 % normal
                hp=errorbar(X,Y,YSTD,'color',Col,'linewidth',obj.p_line_width,...
                    'Marker',obj.p_marker,'MarkerSize',obj.p_marker_size,'LineStyle',obj.p_line_style); hold on;
                %hp=plot(X, Y, obj.p_line_style, 'Color', Col,'linewidth',obj.p_line_width,...
                %'Marker',obj.p_marker,'MarkerSize',obj.p_marker_size);
            elseif Shaded==2 %bar
                hp=bar(X, Y,obj.bar_width,'FaceColor', Col); hold on;
                errorbar(X,Y,YSTD,'.','color',obj.ErrBarColor);
           %     arrayfun(@(x) text(X(x)-.3,Y(x)+YSTD(x)+1,num2str(Y(x),2)),1:length(X));
            elseif Shaded==3 % just plot it
                hp=plot(X, Y, obj.p_line_style, 'Color', Col,'linewidth',obj.p_line_width,...
                    'Marker',obj.p_marker,'MarkerSize',obj.p_marker_size); hold on;
                % now add a merker every N Trials
                data_counter = 0;
                % Plot the data points with shading and markers at specific intervals
                for i = 1:length(X)-1
                    data_counter = data_counter + 1;
                    if data_counter == obj.p_MarkerNpnts
                        % Add a marker
                        plot(X(i), Y(i), 'Color', Col,'Marker',obj.p_marker,'MarkerSize',obj.p_marker_size);
                        data_counter = 0;  % Reset the counter
                    end
                end
            elseif Shaded==4 % plot without line
                hp=plot(X, Y, '.', 'Color', Col,'LineStyle',obj.p_line_style,'linewidth',obj.p_line_width,...
                    'Marker',obj.p_marker,'MarkerSize',obj.p_marker_size,'MarkerFaceColor',Col,'MarkerEdgeColor',Col);
            elseif Shaded==5   % plot with dotted lines             
                hp=plot(X,Y,'LineStyle',obj.p_line_style, 'Color', Col,'linewidth',obj.p_line_width,'Marker',obj.p_marker,'MarkerSize',obj.p_marker_size); hold on;
                plot(X,Y + YSTD,'LineStyle',':','Color', Col,'linewidth',obj.p_line_width,'Marker',obj.p_marker,'MarkerSize',obj.p_marker_size); 
                plot(X,Y - YSTD,'LineStyle',':', 'Color', Col,'linewidth',obj.p_line_width,'Marker',obj.p_marker,'MarkerSize',obj.p_marker_size);                            
            elseif Shaded==6 % plot shaded error bar with a line marker that is repeting every N points
                shadedErrorBar(X,Y,YSTD,'lineProps',{'color',Col,'LineStyle',obj.p_line_style,...
                    'linewidth',obj.p_line_width,'Marker','none'},'patchSaturation',0.2); 
                alpha(obj.patch_alpha);
                hold on;
                hp=plot(X, Y, 'LineStyle',obj.p_line_style, 'Color', Col,'linewidth',obj.p_line_width);                 
                hp.Visible='off';
                % now add a merker every N Trials
                data_counter = 0;
                % Plot the data points with shading and markers at specific intervals
                for i = 1:length(X)-1
                    data_counter = data_counter + 1;
                    if data_counter == obj.p_MarkerNpnts
                        % Add a marker
                        plot(X(i), Y(i), 'Color', Col,'Marker',obj.p_marker,'MarkerSize',obj.p_marker_size);
                        data_counter = 0;  % Reset the counter
                    end
                end
            end
            xlabel(Xlbl);ylabel(Ylbl);
            grid(obj.AxesGrid);
            %obj.AppendFigTitles(gca,[Title; {[' n=' num2str(n)]}]);
            obj.AppendLegend(obj.LegendTxt,[],hp,ActObj,ActObjLegends); % add legend to current legends
            obj.AddTitle(Title,n);
           
            % if X is s time axis then autmatically set the axis limits
            if sum(contains(Xlbl,'Time','IgnoreCase',1))
                obj.SetBothTimeAxisTicks('x','auto',AnalysisOpts.PaperSpec.TimeAxisSteps,AnalysisOpts.PaperSpec.TimeAxisSteps,'auto',[]);
            end
            obj.FormatAxes(gca);

        end
        function plotDataWithShadingAndMarkers(obj,x, y,YSTD,Col, N,varargin)

            global AnalysisOpts
            obj=obj.ParseParams(varargin);

            % Check if N is not provided as an input, use a default value
            if isempty(N) 
                N = 5;
            end
            hold on
            % Plot the data points
            plot(x, y,Col,'linewidth',obj.p_line_width,...
                    'LineStyle',obj.p_line_style);

            % Initialize arrays for shaded patches
            xPatch = [];
            yPatch = [];

            % Initialize a counter to keep track of the data point index
            data_counter = 0;

            % Define the shading color and transparency
            shading_color = Col;  % Blue
            alpha = 0.2;  % Transparency for shading

            % Plot the data points with shading and markers at specific intervals
            for i = 1:length(x)-1
                data_counter = data_counter + 1;

                if data_counter == N
                    % Add a marker
                    plot(x(i), y(i), 'Marker',obj.p_marker,'MarkerSize',obj.p_marker_size);
                    data_counter = 0;  % Reset the counter

                end

                    % Add the shaded area for the past N data points
                    xPatch = [xPatch, x(i), (x(i+1))];
                    yPatch = [yPatch, y(i) + YSTD(i), (y(i+1) - YSTD(i+1))];

                    % Create the shaded area
                    patch(xPatch, yPatch, shading_color, 'FaceAlpha', alpha, 'EdgeColor', 'none');
            end          

        end

        function AddTitle(obj,Title,n) % adds title
            Title=obj.ManData.RepDashSpace(Title);
            if obj.include_n & ~isempty(n)
                T=title([Title;{[' n=' num2str(n)]}]);
            else
                T=title(Title);
            end
            T.FontSize=obj.TitleFontSize;
        end
        function AddSignificanceStar(obj,X,hsig,Col,Sp,varargin) % adds significance star to a plot 
            obj=obj.ParseParams(varargin) ; %%Process optional inputs
            if min(size(X))>1;return;end % more than one pvalue per raw and can't be processed
            if ~isempty(Sp)
               subplot(Sp);
            end
            v=axis;ymax=v(4);
            text(X(hsig>0),ymax*ones(1,sum(hsig>0)),'*','Color',Col,'FontSize',obj.SigStar_fontsize)
            ylim([v(3) ymax+0.05*ymax])
        end
         function AddDetailedSignificanceStar(obj,X,pval,Col,Sp,varargin) % adds significance star to a plot 
            obj=obj.ParseParams(varargin) ; %%Process optional inputs
            if min(size(X))>1;return;end % more than one pvalue per raw and can't be processed
            if ~isempty(Sp)
               subplot(Sp);
            end
            v=axis;ymax=v(4);ymin=v(3);D=ymax-ymin;
            Symb=pvalueStar(pval);
            text(X,ymax+0.05*D,Symb,'Color',Col,'FontSize',obj.SigStar_fontsize)
            text(X,ymax+0.1*D,num2str(pval,4),'Color',Col,'FontSize',obj.SigStar_fontsize)
         %   ylim([v(3) ymax+D])
        end
        function [TimeInd,StrTim,EndTim]=LimitTimeAxis(obj,X,xlbl,varargin)
            global AnalysisOpts
            obj=obj.ParseParams(varargin) ; %%Process optional inputs
            % if this is not a time a time axis then don't use it
            if ~obj.IsthisAxisTime;TimeInd=logical(ones(size(X)));StrTim=[];EndTim=[];return;end
           
            if ~isempty(AnalysisOpts.ThisTimeAxisStart) && ~isempty(AnalysisOpts.ThisTimeAxisEnd)
                StrTim=AnalysisOpts.ThisTimeAxisStart;
                EndTim=AnalysisOpts.ThisTimeAxisEnd;
            else
                StrTim=AnalysisOpts.PaperSpec.(['StrTime_' AnalysisOpts.SpkCntStartFieldName]);
                EndTim=AnalysisOpts.PaperSpec.(['EndTime_' AnalysisOpts.SpkCntStartFieldName]);
            end
            % if we are limiting time axis then check if the label has time in it then limit
            if AnalysisOpts.PaperSpec.LimitTimeAxis & sum(contains(xlbl,'Time','IgnoreCase',1))                
                TimeInd=X>=StrTim & X<=EndTim;               
            else
                TimeInd=logical(ones(size(X)));
            end
        end
        function SetAxisLimits(~,AxisType,AxisLimits,Dim) % set limit for c,x,y,x axes
            % Dim is the raw where we use the Axis Limits matrix
            if ~strcmp(AxisLimits,'auto')
                AxisLimits=AxisLimits(Dim,:);
            end
            switch AxisType
                case 'c'
                    caxis(AxisLimits);
                case {'x','y','z'}
                    eval(sprintf('%slim(AxisLimits);',AxisType));
            end                
        end
        function SetTimeAxisLims(obj) % limits the current time axis
            global AnalysisOpts
           
            if ~AnalysisOpts.PaperSpec.LimitTimeAxis;return;end
            if ~isempty(AnalysisOpts.ThisTimeAxisStart) && ~isempty(AnalysisOpts.ThisTimeAxisEnd)
                StrTim=AnalysisOpts.ThisTimeAxisStart;
                EndTim=AnalysisOpts.ThisTimeAxisEnd;
            else
                StrTim=AnalysisOpts.PaperSpec.(['StrTime_' AnalysisOpts.SpkCntStartFieldName]);
                EndTim=AnalysisOpts.PaperSpec.(['EndTime_' AnalysisOpts.SpkCntStartFieldName]);
            end
            
            xlim(gca,[StrTim EndTim])            
        end
        function [hp]=PolarPlot(obj,theta,rho,thetaTicks,Col,Title,varargin)
            if ~isempty(varargin) & numel(varargin)==1
                LineProp=varargin{1};
            else
                LineProp='-';                
                obj=obj.ParseParams(varargin) ; %%Process optional inputs
            end
            n=size(rho,1); % number of samples
            
            Col=obj.getSingleColor(Col);
            if ~isempty(obj.ThisSubplot);subplot(obj.ThisSubplot(1),obj.ThisSubplot(2),obj.ThisSubplot(3));end           
           
            
            % polar plot of data 
            hp=polarplot(theta,rho,'LineStyle',LineProp,'Color', Col,'linewidth',obj.p_line_width,'Marker',obj.polar_marker,'MarkerFaceColor','auto');
            pax=gca;
            pax.ThetaAxisUnits = 'radians';
            thetaticks(theta);
            thetaticklabels(thetaTicks);
            pax.ThetaColor = obj.polar_ThetaColor;
            pax.RColor =  obj.polar_RColor;
            title([Title ]);%' n=' num2str(n)]);
            obj.FormatPolarAxes(gca);
        end
        function [hp]=Plot(obj,X,Y,Col,Xlbl,Ylbl,Title,Sp,varargin) % plots whatever
            obj=obj.ParseParams(varargin) ; %%Process optional inputs
            Col=obj.getSingleColor(Col); 
            if ~isempty(obj.ThisSubplot);subplot(obj.ThisSubplot(1),obj.ThisSubplot(2),obj.ThisSubplot(3));end
            if ~isempty(Sp)
               subplot(Sp);
            end
            
            hp=plot(X,Y,obj.p_line_style,'linewidth',obj.p_line_width,'MarkerSize',obj.p_marker_size,...
                'color',Col);
            xlabel(Xlbl);ylabel(Ylbl);title(Title);
            % format all axis
            obj.FormatAxes(gca);
        end
        function [hp]=ScatterPlot(obj,X,Y,Col,Xlbl,Ylbl,Title,Sp,RemoveOutLiersFlag,varargin) % scaterplots whatever
            obj=obj.ParseParams(varargin) ; %%Process optional inputs
            Col=obj.getSingleColor(Col); 
            if ~isempty(obj.ThisSubplot);subplot(obj.ThisSubplot(1),obj.ThisSubplot(2),obj.ThisSubplot(3));end
            if ~isempty(Sp)
               subplot(Sp);
            end
            if RemoveOutLiersFlag % remove outliers from  the data
                [X,Y]=obj.ManData.removeoutliers(X,Y,[]);
            end
            hp=scatter(X,Y,obj.p_marker_size,'filled', 'MarkerFaceColor',Col);
            xlabel(Xlbl);ylabel(Ylbl);title(Title);
            % format all axis
            obj.FormatAxes(gca);
        end
        
        function [hp]=BarPlot(obj,X,Y,Col,Xlbl,Ylbl,Title,varargin) % plots whatever
            obj=obj.ParseParams(varargin) ; %%Process optional inputs
            Col=obj.getSingleColor(Col);
            if ~isempty(obj.ThisSubplot);subplot(obj.ThisSubplot(1),obj.ThisSubplot(2),obj.ThisSubplot(3));end
                                   
            [ActObj,ActObjLegends]=obj.GetActiveLegends([]);
                 
            hp=bar(X,Y,'FaceColor',Col,'EdgeColor',obj.bar_edgecol,'FaceAlpha',obj.bar_face_alpha,...
                'EdgeAlpha',obj.bar_edge_alpha,'BarWidth',obj.bar_width);        
            xlabel(Xlbl);ylabel(Ylbl);
            obj.AddTitle(Title,[]); 
             
            obj.AppendLegend(obj.LegendTxt,[],hp,ActObj,ActObjLegends); % add legend to current legends
            % format all axis
            obj.FormatAxes(gca);
        end
        function hp=HistogramPlot(obj,X,edges,Col,Xlbl,Ylbl,Title,varargin) % plots histogram 
            obj=obj.ParseParams(varargin) ; %%Process optional inputs
            Col=obj.getSingleColor(Col);
            if ~isempty(obj.ThisSubplot);subplot(obj.ThisSubplot(1),obj.ThisSubplot(2),obj.ThisSubplot(3));end
            
            if isempty(edges)
                [N,edges] = histcounts(X);
            else
                [N,edges] = histcounts(X,edges);
            end
            if obj.histplot_calperc % are we showing proportion 
                N=N/sum(N);
            end
             hp=obj.BarPlot(movmean(edges,2,'Endpoints','discard'),N,Col,Xlbl,Ylbl,Title);
            
            % add the mean value there 
            hold on
            MeanN=nanmean(X);
            v=axis;
            plot([MeanN MeanN],[v(3) v(4)],'--r','linewidth',obj.p_line_width);
             obj.FormatAxes(gca);
        end
        function hp=Text(obj,X,Y,txt,Col,varargin)
            obj=obj.ParseParams(varargin) ; %%Process optional inputs
            Col=obj.getSingleColor(Col);
            if ~isempty(obj.ThisSubplot);subplot(obj.ThisSubplot(1),obj.ThisSubplot(2),obj.ThisSubplot(3));end
            
            hp=text(X,Y,txt,'color',Col,'FontSize',obj.font_size,'FontName',obj.font_name,...
                'FontWeight',obj.font_weight,'FontUnits',obj.unitss,'FontAngle',obj.font_angle);
            
        end
        %% table functions
        function handles = TableFigure(obj,T,handles)
            
            %if no figure handle is passed, then plots new figure
            if nargin ==2
                handles = figure;
            end
            % Get the table in string form.
            format shortG
            TString = evalc('disp(T)');
            % Use TeX Markup for bold formatting and underscores.
            TString = strrep(TString,'<strong>','\bf');
            TString = strrep(TString,'</strong>','\rm');
            TString = strrep(TString,'_','\_');
            % Get a fixed-width font.
            FixedWidth = get(0,'FixedWidthFontName');
            % Output the table using the annotation command.
            annotation(handles,'Textbox','String',TString,'Interpreter','Tex',...
                'FontName',FixedWidth,'Units','Normalized','Position',[0 0 1 1])
            
            %Note: alternative is to use
            % uitable('Data',T{:,:},'ColumnName',T.Properties.VariableNames,...
            %     'RowName',T.Properties.RowNames,'Units', 'Normalized', 'Position',[0, 0, 1, 1],...
            %     'FontWeight','bold','FontSize',12);
            
        end
        function done=plot_significance_level(obj,X,P,Time,AxisLims,c,thr,MinYAxisLim,varargin) % plot significance level 
            % @x  axis indices
            % @ Time x axis x ticks
            % @p p value for each time point 
            % @AxisLims current axis limit; if you give one number then it
            % puts the levels there.
            % @c color of the significance level 
            % @thr p value thresholds. if empty it would be [0.05 0.01 0.001]
            % @MinAxisLim minimum of axis limits
            global AnalysisOpts
            
            obj=obj.ParseParams(varargin) ; %%Process optional inputs
           
            if ~isempty(obj.Sp)
                subplot(obj.Sp);
            end
            if strcmp(AxisLims,'auto');axis tight;v=axis;AxisLims=v(3:4);end
            if length(AxisLims)>=2
                v=axis;
                ThisAxisLimit=max(v(4),AxisLims(2));%ceil(v(4) * 10) / 10 - 0.05;%
                if length(AxisLims)==3
                    DeltaYlim=AxisLims(3);
                else
                    % change the ylimits so that it matches
                    DeltaYlim=abs(AxisLims(2)-AxisLims(1))/30;
                end
                a=ThisAxisLimit+DeltaYlim;
                ylim([v(3) ThisAxisLimit+2*DeltaYlim]);                           
            elseif length(AxisLims)==1
                a=AxisLims;
            end
            if isempty(X) & isempty(P);done=1;return;end

            if isempty(thr)
               thr=[0.05 0.01 0.001];
            end
            if isempty(c)
                c=obj.getSingleColor(1);
            end         
%             if obj.SmoothTimeAxis & obj.IsthisAxisTime% smooth Time axis                
%                 Time=obj.ManData.SmoothData(Time,obj.WidthSmoothing,'SmoothingMethod',obj.SmoothingMethod);
%             end
            % are we limiting the time axis
            if AnalysisOpts.PaperSpec.LimitTimeAxis & obj.IsthisAxisTime
                [TimeInd,StrTim,EndTim]=obj.LimitTimeAxis(Time,'Time');
                xlim([StrTim EndTim]);
                if ~iscell(P) & length(Time)==length(P) % then we have the pvalues themselves 
                    % now apply correction to this 
                    Pcorr=obj.ManData.PerformStattestMultipleComparisionCorrection(P(TimeInd));
                  %  figure;plot(Time,P);hold on;
                    P(TimeInd)=Pcorr;
                  %  plot(Time,P);
                  %  legend('UnCorrected pval','Corrected pval');
                    [X,P]=obj.ManData.DifferentiateSigClusters2(P);                    
                elseif iscell(P) & isempty(X) % then we are dealing with  trl shuffle
                    nSigLvl=length(AnalysisOpts.Modified_MannKendall_significance_value_ac);
                    PindTot=nan*ones(nSigLvl,sum(TimeInd));
                    % get different values
                    for i=1:nSigLvl
                        Pind{i}=obj.ManData.PerformStattestMultipleComparisionCorrection(cellfun(@(x) x(i),P(TimeInd)));
                        PindTot(i,Pind{i}<0.05)=AnalysisOpts.Modified_MannKendall_significance_value_ac(i)-0.0001;
                        PindTot(i,Pind{i}>=0.05)=1;
                    end
                    PindMin=min(PindTot,[],1);
                    PindF=nan*ones(1,length(P));
                    PindF(TimeInd)=PindMin;
                    [X,P]=obj.ManData.DifferentiateSigClusters2(PindF); 
                end
                TimeX=cellfun(@(x) Time(x),X,'UniformOutput',0);
                [~,~,XInc]=cellfun(@(x) obj.ManData.RemoveEntryFromVec(x,x(x<StrTim | x>EndTim)),TimeX,'UniformOutput',0);
                X=arrayfun(@(x) X{x}(XInc{x}),1:length(XInc),'UniformOutput',0);
                P=arrayfun(@(x) P{x}(XInc{x}),1:length(XInc),'UniformOutput',0);
            end
            
            if ~iscell(X);X={X};P={P};end
            
            % before adding get the active legends in case we want to add
            % legends in the future 
            [ActObj,ActObjLegends]=obj.GetActiveLegends([]);           
            hold on
            for nc=1:length(X) % loop on clusters
                x=X{nc};p=P{nc};
                if ~isempty(x)
                    for b=1:length(thr)
                        this_thr=thr(b);
                        
                        for i=1:length(p)-1
                            if i==1 & AnalysisOpts.ShowStatPvalinPlot; text(Time(x(i)),a-0.05*a,num2str(p(i),4));end
                            if p(i)<this_thr && p(i+1)<this_thr
                                plot(Time(x(i):x(i+1)),a*ones(1,length(x(i):x(i+1))),'-','Color',c,'LineWidth',obj.sig_line_width(b));
                            end
                            if i>1
                                if p(i-1)>=this_thr && p(i)<this_thr && p(i+1)>=this_thr
                                    plot(Time(x(i)),a,'.','Color',c,'LineWidth',obj.sig_line_width(b));
                                end
                            elseif i==1
                                if  p(i)<this_thr && p(i+1)>=this_thr
                                    plot(Time(x(i)),a,'.','Color',c,'LineWidth',obj.sig_line_width(b));
                                end
                            end
                        end
                    end
                end
            end
            if exist('MinYAxisLim','var') % if we want to put limit of min Y axis limit
                if ~isempty(MinYAxisLim)
                    v=axis;
                    ylim([MinYAxisLim v(4)]);
                    
                end
            end 
            
           
            obj.AppendLegend([],[],[],ActObj,ActObjLegends) ; % keep whatever we have here  
            done=1;% we have applied this significance
        end
        function plot_significance_levelNoCluster(obj)
            % Example data (replace with your actual data)
            time = 1:100;
            p_values = rand(1, 100); % Replace with your p-values

            % Create a figure and plot your data
            figure;
            plot(time, p_values);

            % Set the significance levels
            significance_levels = [0.001, 0.01, 0.05];

            % Loop through each significance level and add line segments for islands of significance
            for significance_level = significance_levels
                % Find time points where p-values are below the significance level
                significant_indices = find(p_values < significance_level);

                if ~isempty(significant_indices)
                    island_start = significant_indices(1);

                    for i = 2:length(significant_indices)
                        if significant_indices(i) ~= significant_indices(i - 1) + 1
                            % Add a line segment for the island of significance
                            island_end = significant_indices(i - 1);
                            line_width = 2; % You can adjust the line width
                            line_color = 'r'; % You can choose a color
                            y = max(p_values) + 0.05; % Adjust the vertical position of the lines
                            line([time(island_start), time(island_end)], [y, y], ...
                                'Color', line_color, 'LineWidth', line_width);

                            % Update the island_start
                            island_start = significant_indices(i);
                        end
                    end

                    % Add a line segment for the last island of significance
                    island_end = significant_indices(end);
                    line_width = 2; % You can adjust the line width
                    line_color = 'r'; % You can choose a color
                    y = max(p_values) + 0.05; % Adjust the vertical position of the lines
                    line([time(island_start), time(island_end)], [y, y], ...
                        'Color', line_color, 'LineWidth', line_width);
                end
            end

            % Adjust axis limits and labels if needed
            xlabel('Time');
            ylabel('P-values');

            % Customize your plot as needed


        end
        function hp=PlotVerticalLine(obj,t,Sp,Col,varargin) % plots a vertical line on the current plot
            % t is time point we want the bar to be ploted
            obj=obj.ParseParams(varargin) ; %%Process optional inputs

            if ~isempty(obj.ThisSubplot);subplot(obj.ThisSubplot(1),obj.ThisSubplot(2),obj.ThisSubplot(3));end
            if ~isempty(Sp)
               subplot(Sp);
            end
            [ActObj,ActObjLegends]=obj.GetActiveLegends([]);
            Col=obj.getSingleColor(Col);
            v=axis;
            hp=plot([t t],[v(3) v(4)],obj.p_line_style,'linewidth',obj.OriginLine_width,...
                'MarkerSize',obj.p_marker_size,'color',Col);
            obj.AppendLegend(obj.LegendTxt,[],hp,ActObj,ActObjLegends); % add legend to current legends
        end
        
        function hp=PlotHorizontalLine(obj,t,Sp,Col,varargin)
            global AnalysisOpts
            obj=obj.ParseParams(varargin) ; %%Process optional inputs

            if ~isempty(obj.ThisSubplot);subplot(obj.ThisSubplot(1),obj.ThisSubplot(2),obj.ThisSubplot(3));end
            if ~isempty(Sp)
               subplot(Sp);
            end
            Col=obj.getSingleColor(Col);
            v=axis;
            [ActObj,ActObjLegends]=obj.GetActiveLegends([]);

            hp=plot([v(1) v(2)],[t t],obj.p_line_style,'linewidth',obj.OriginLine_width,'MarkerSize',obj.p_marker_size,...
                'color',Col);       
            obj.AppendLegend(obj.LegendTxt,[],hp,ActObj,ActObjLegends); % add legend to current legends

        end
        function [h1,TimeIndX,TimeIndY,X,Y,Z]=Image(obj,X,Y,Z,Xlbl,Ylbl,ColBarlbl,Title,Sp,varargin) % plots an image 
            obj=obj.ParseParams(varargin) ; %%Process optional inputs
            
            if ~isempty(Sp);axes(Sp);end
            % smoothdata if we need to 
            if ~isempty(obj.SmoothingMethod)    
                % if both axis are time then use convolution to smooth
                if sum(contains(Xlbl,'Time','IgnoreCase',true)) &  sum(contains(Ylbl,'Time','IgnoreCase',true))
                    % Z = convn(Z, ones(obj.WidthSmoothing)/(obj.WidthSmoothing^2), 'same');
                    % kernel = fspecial('average', [obj.WidthSmoothing obj.WidthSmoothing]);
                    % Z=imfilter(Z,kernel,'replicate');
                    % Z=smoothdata2(Z,"movmean",obj.WidthSmoothing);
                    Z=obj.ManData.SmoothData(Z,obj.WidthSmoothing,'SmoothingMethod',obj.SmoothingMethod,'DimSmoothing',1);
                    Z=obj.ManData.SmoothData(Z,obj.WidthSmoothing,'SmoothingMethod',obj.SmoothingMethod,'DimSmoothing',2);
          
                else %otherwise first smooth the rows and then columns
                    Z=obj.ManData.SmoothData(Z,obj.WidthSmoothingDim2,'SmoothingMethod',obj.SmoothingMethod,'DimSmoothing',1);
                    Z=obj.ManData.SmoothData(Z,obj.WidthSmoothing,'SmoothingMethod',obj.SmoothingMethod,'DimSmoothing',2);
                end
                if obj.SmoothTimeAxis & obj.IsthisAxisTime% if we are smoothing the time axis as well
                    % smooth time axis too, we assume that columns represent time
                    if sum(contains(Xlbl,'Time','IgnoreCase',true))
                        X=obj.ManData.SmoothData(X,obj.WidthSmoothing,'SmoothingMethod',obj.SmoothingMethod,'DimSmoothing',2);
                        warning('X axis might not be time but it has been smoothed');
                    end
                    
                    if sum(contains(Ylbl,'Time','IgnoreCase',true));
                        warning('Y axis contains time and has been smoothed');
                        Y=obj.ManData.SmoothData(Y,obj.WidthSmoothing,'SmoothingMethod',obj.SmoothingMethod,'DimSmoothing',2);
                    end
                end
            end
            % change time axis if we need to 
            TimeIndX=obj.LimitTimeAxis(X,Xlbl); 
            TimeIndY=obj.LimitTimeAxis(Y,Ylbl);            
            X=X(TimeIndX);Y=Y(TimeIndY);Z=Z(TimeIndY,TimeIndX);
            
            if strcmp(obj.imageplotfunc,'surf')
                args = {X,Y,Z};
                h1=surf(args{:},'edgecolor','none');
                view(0,90);
              %  shading interp;
            elseif strcmp(obj.imageplotfunc,'imagesc')
                h1=imagesc(X,Y,Z);
                set(gca,'YDir','normal');
            elseif strcmp(obj.imageplotfunc,'pcolor')
                 h1=pcolor(X,Y,Z);
                 h1.EdgeColor='none';
                 h1.FaceColor = 'interp';
                 set(gca,'YDir','normal')
            end
            
            axis tight;
            colormap(obj.image_colormap);
            h = colorbar;
            caxis(obj.caxis_limits);
            h.Label.String = ColBarlbl;
            xlabel(Xlbl); 
            ylabel(Ylbl);
            obj.AddTitle(Title,[]);
                        
            if ~isempty(obj.Xticks) | ~isempty(obj.Yticks)
                xticks(obj.Xticks);
                xticklabels(obj.XticksLabels);
                yticks(obj.Yticks);
                yticklabels(obj.YticksLabels);               
            end
            % add origin line 
            obj.addOriginLines2Image(obj.OriginLine,Z)
            % format all axis
            obj.FormatAxes(gca);
        end
        
        function addOriginLines2Image(obj,OriginLine,Z)
            hold on
            v=axis;
            switch OriginLine
                case 1 % if we are plotting a line at the [0 0]
                    plot3([v(1) v(2)],[0 0],[max(Z(:)) max(Z(:))],'Color',[1 1 1],'LineWidth',obj.OriginLine_width);
                    plot3([0 0],[v(3) v(4)],[max(Z(:)) max(Z(:))],'Color',[1 1 1],'LineWidth',obj.OriginLine_width);
                case 2 % only show on the x-axis
                    plot3([0 0],[v(3) v(4)],[max(Z(:)) max(Z(:))],'Color',[1 1 1],'LineWidth',obj.OriginLine_width);
                case 3 % only show on the y-axis
                    
                    plot3([v(1) v(2)],[0 0],[max(Z(:)) max(Z(:))],'Color',[1 1 1],'LineWidth',obj.OriginLine_width);
                case 4 % line on x=y
                    plot3([v(1) v(2)],[v(3) v(4)],[max(Z(:)) max(Z(:))],'Color',[1 1 1],'LineWidth',obj.OriginLine_width);
                case 5 % line on x=y  and line at [0 0]
                    plot3([v(1) v(2)],[0 0],[max(Z(:)) max(Z(:))],'Color',[1 1 1],'LineWidth',obj.OriginLine_width);
                    plot3([0 0],[v(3) v(4)],[max(Z(:)) max(Z(:))],'Color',[1 1 1],'LineWidth',obj.OriginLine_width);
                    plot3([v(1) v(2)],[v(3) v(4)],[max(Z(:)) max(Z(:))],'Color',[1 1 1],'LineWidth',obj.OriginLine_width);
            end
        end
        function AppendFigTitles(obj,h,Title,varargin) % adds a title to the correct title of the paper if there is already one
             obj=obj.ParseParams(varargin) ; %%Process optional inputs

            if isprop(h.Children,'Title') & obj.AppendTitles % if we aleady have a title add this one in the second line 
                CurrTitle=get(h,'Title');
                set(h,'Title',[CurrTitle.String;{Title}]);
            else
                title(Title);
            end
        end
        function AppendLegend(obj,LegTxt,h,hAdd,ActObj,ActObjLegends)
            hold on
            % h is the axis object 
            % hAdd is the object that needs to be added to the current plot
            % ActObj previous active objects that had legend
            % Legend of previous active objects
            if isempty(LegTxt) & isempty(ActObj)                
                return; 
            elseif isempty(LegTxt) & ~isempty(ActObj) % we want to keep current legends without adding anymore
                obj.FlushLegends([ActObj]);
                return
            end
            if isempty(h);h=gca;end
            
            % get current values of legend
            if isempty(hAdd)
                if isempty(h.Legend) % generate a new legend
                    legend(LegTxt,'Location',obj.LegendLoc,'FontSize',obj.LegendFontSiz);
                    legend('boxoff');
                else
                    h.Legend.String{end}=LegTxt;
                end   
            else
                legend(h,[ActObj; hAdd],[ActObjLegends LegTxt],'Location',obj.LegendLoc,'FontSize',...
                    obj.LegendFontSiz); % adds the curent legend to the prevous active objects 
                legend('boxoff');
                obj.FlushLegends([ActObj;hAdd]) % empty non used legends 
            end
        end
        function [ActObj,ActObjLegends]=GetActiveLegends(~,h)% gets currect active objects with legends and the value of their legends
            if isempty(h);h=gca;end
            
            Child=get(gca,'Children');
            ActiveObjsInds=arrayfun(@(x) ~isempty(Child(x).DisplayName),1:length(Child));
            ActObj=Child(ActiveObjsInds);
            ActObjLegends=arrayfun(@(x) Child(x).DisplayName,find(ActiveObjsInds),'UniformOutput',0);
        end
        function FlushLegends(obj,ActObjs) % flushes the legends that are not used to empty
            % Function to organize legends and hide titles for specified active objects
            % Inputs:
            % - obj: An instance of a class (or struct) containing settings
            % - ActObjs: An array of handles to specific plot elements to work with
            
            % Get all child objects (plot elements) of the current axes
            Child = get(gca, 'Children');
            
            % Find the index of these active objects in the Child array
            % (Assuming ActObjs is an array of handles to specific plot elements you want to work with)
            IndActObj = cell2mat(arrayfun(@(x) find(Child == ActObjs(x)), 1:length(ActObjs), 'UniformOutput', 0));
            
            % Find indices of non-active objects
            NonActiveObj = setdiff(1:length(Child), IndActObj);
            
            % Empty the display names (titles) from non-active objects
            for i = NonActiveObj
                Child(i).DisplayName = '';
            end
            
            % Get display names (titles) of active objects for creating legends
            ActObjLegends = arrayfun(@(x) Child(x).DisplayName, IndActObj, 'UniformOutput', 0);
            
            % Get the current axes handle
            h = gca;
            
            % Create a legend for active objects with the extracted legends and other properties
            legend(Child(IndActObj), ActObjLegends, 'Location', h.Legend.Location, 'FontSize', obj.LegendFontSiz);
            
            % Turn off legend's bounding box
            legend('boxoff');
        end
        %% movie making functions 
        function MakeMovieFromFrames(obj,mvFrame,fps,saveMovieStr)
            % create move from the frames
            % inputs
            % MvFrame       Frames captured with     mvFrame(i) = getframe(gcf);
            %fps            Frame per second
            %saveMovieStr   Save Path
            % use     mvFrame(cur_win) = getframe(gcf);
            
            saveMovieStr = [saveMovieStr, '.avi'];
            
            writerObj = VideoWriter(saveMovieStr, 'Motion JPEG AVI');
            writerObj.FrameRate = fps;
            open(writerObj);
            fprintf('writing movie...\n')
            mvFrame = mvFrame(~cellfun(@isempty,{mvFrame.cdata})); %remove dropped frames
            try
                for i =1:length(mvFrame)
                    if isempty(mvFrame(i).cdata)
                        mvFrame(i) = [];
                    end
                end
            catch
            end
            writeVideo(writerObj, mvFrame);
            close(writerObj);
        end
        function SetAxisTicks(obj,AxisName,Lims,Steps)
            if strcmp(Lims,'auto')
                Lims=axis;
                LimsX=Lims(1:2);LimsY=Lims(3:4);
            else
                LimsX=Lims;LimsY=Lims;
            end
            TicksX=LimsX(1):Steps:LimsX(2);
            TicksY=LimsY(1):Steps:LimsY(2);
            switch AxisName 
                case 'x'
                    xticks(TicksX);
                case 'y'
                    yticks(TicksY);
            end            
        end
        function SetBothTimeAxisTicks(obj,AxisName,Lims,StepsX,StepsY,TicksXLabels,TicksYLabels)
            global AnalysisOpts
            [TimeInd,StrTim,EndTim]=obj.LimitTimeAxis(AnalysisOpts.Time,'Time','IsthisAxisTime',1);
            if strcmp(Lims,'auto')
                 LimsX=[StrTim EndTim];LimsY=[StrTim EndTim];
            else
                LimsX=Lims(1,:);LimsY=Lims(2,:);
            end
            TicksX=LimsX(1):StepsX:LimsX(2);
            TicksY=LimsY(1):StepsY:LimsY(2);
            switch AxisName 
                case 'x'
                    xticks(TicksX);
                    if strcmp(TicksXLabels,'auto')
                        xticklabels(obj.ManData.CovertDouble2CellStr(TicksX));
                    elseif ~isempty(TicksXLabels)
                        xticklabels(TicksXLabels);
                    end
                case 'y'
                    yticks(TicksY); 
                    if strcmp(TicksYLabels,'auto')
                        yticklabels(obj.ManData.CovertDouble2CellStr(TicksY));
                    elseif ~isempty(TicksYLabels)
                        yticklabels(TicksYLabels);
                    end
            end            
        end
        function shiftedColormap=ShiftColorMap(obj,CurrentColorMap,caxisLimits,Center)
            %In this code, the colormap values are shifted to center the colormap on Center=0.5 
            %while maintaining the specified color axis limits of 0.3 to 0.7. 
            %The centerIndex is calculated based on the desired center value (0.5) and the color axis limits.
            %Then, a portion of the original colormap is placed around the center index to create the shifted colormap. 
            %Finally, the shifted colormap is visualized using the imagesc function.           
            %Remember to replace jet(256) with your actual colormap and adjust the visualization settings as needed.

            % Define the colormap
            colormapValues = CurrentColorMap;%jet(256);  % Replace with your actual colormap
            
            % Define the color axis limits
            %caxisLimits = [0.3, 0.7];
            centerIndex = round((Center *(size(colormapValues, 1) - 1)) + 1);
            UpperNumIndeces=round(((Center-caxisLimits(1)) *(size(colormapValues, 1) - 1)) + 1);
            LowerNumberIndeces=round(((caxisLimits(2)-Center) *(size(colormapValues, 1) - 1)) + 1);
           
            % Calculate the range of indices to shift
            shiftRange = centerIndex - UpperNumIndeces:centerIndex + LowerNumberIndeces;
            
            % Create a new colormap that is centered on 0.5
            shiftedColormap = colormapValues(shiftRange, :);
%             % Create a figure to visualize the colormap
%             figure;
%             imagesc([0.3, 0.7], [1, 256], shiftedColormap);
%             colorbar;
            colormap(shiftedColormap);
            % Set the color axis limits
            caxis(caxisLimits);
            
        end
       
        function CopyFigures(obj,Ax1,Ax2,YLIM,Index)
            h=figure;ax_combined=gca;
            copyobj(get(Ax1,'children'),ax_combined);
            copyobj(get(Ax2,'children'),ax_combined);
            axis tight;
            ylim(YLIM);
            currentDateTime = clock;
            
            % Convert numbers to characters
            dateString = sprintf('%04d%02d%02d', currentDateTime(1:3));
            if isempty(Index);Index=randi(100);end
            [~,~,FigFileName]=obj.ManData.GetFileName(['Classifier'],['CombinedFig_' dateString '_' num2str(Index)],'SaveInResults',1,'WantedDate','ALL');
            obj.SaveFigSeries([],FigFileName,[h],'SaveEachFrame',1,'SaveEachSubplots',1);
            
        end
    end %methods
    
end %class structure
