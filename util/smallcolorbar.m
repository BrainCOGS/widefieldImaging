function ch = smallcolorbar(axs,location,axiscl)

% ch = smallcolorbar(axs,location,axiscl)
% plots a small colorbar that doesn't squish the axis
% axs is axis handle
% location is 'eastoutside' (default), 'southoutside' or 'northoutside'
% axisCl is color of box around colorbar and numbers (default 'k')

if nargin < 1; axs      = gca;           end
if nargin < 2; location = 'eastoutside'; end
if nargin < 3; axiscl   = [0 0 0];       end

axspos = get(axs, 'position');
ratio  = [.08 .33];

switch location
  case 'eastoutside'
    x_size = axspos(3)*ratio(1);
    y_size = axspos(4)*ratio(2);
    pos    = [sum(axspos([1 3])) + x_size*.2, axspos(2), x_size, y_size];
    
  case 'southoutside'
    x_size = axspos(3)*ratio(2);
    y_size = axspos(4)*ratio(1);
    pos    = [axspos(1), axspos(2) - y_size*1.2, x_size, y_size];
  
  case 'northoutside'
    x_size = axspos(3)*ratio(2);
    y_size = axspos(4)*ratio(1);
    pos    = [axspos(1), axspos(2) + axspos(4) + y_size*.1, x_size, y_size];
    
  otherwise
    error('position not recognized')
end

ch = colorbar('location',location,'position',pos,'Color',axiscl);