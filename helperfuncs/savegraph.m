function savegraph(h1,x_text,y_text,varargin)
% SAVEGRAPH: Add formatted text to graph with optional save
%   h1:     plot object
%   x_text:         Text for x axis
%   y_text:         Text for y axis
%   varargin{1}:    Argument passed to legend() function
%   varargin{2}:    File is saved on this name (with location) as jpg

xlhand = get(gca,'xlabel');
set(xlhand,'string',x_text,'fontsize',20);
yrhand = get(gca,'ylabel');
set(yrhand,'string',y_text,'fontsize',20);

if nargin > 3
    h_legend=legend(h1, varargin{1});
    set(h_legend,'FontSize',14);
    if nargin > 4
        saveas(h1,varargin{2},'jpg')
    end
end

end