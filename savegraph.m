function savegraph(h1,x_text,y_text,legend_text,name)

xlhand = get(gca,'xlabel');
set(xlhand,'string',x_text,'fontsize',20);
yrhand = get(gca,'ylabel');
set(yrhand,'string',y_text,'fontsize',20);

% if ~isempty(legend_text)
%     h_legend=legend(legend_text);
%     set(h_legend,'FontSize',14);
% end

saveas(h1,name,'jpg')

end