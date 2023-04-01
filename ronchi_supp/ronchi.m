## Copyright (C) 2020 Khosrow Hassani
## University of Tehran, Iran.
## This program is free software: you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
## 
## This program is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see
## <https://www.gnu.org/licenses/>.

## -*- texinfo -*- 
## @deftypefn {} {@var{g} =} ronchi (@var{nx},@var{ny},@var{d},@va@var{phi},@var{offset},@var{neg}) 
##
## @seealso{ronchigram}
## @end deftypefn

## Author: Khosrow Hassani <hassanikh@ut.ac.ir>
## Created: 2020-11-04

function g=ronchi(nx,ny,d,phi,xoff,neg)
  if nargin<6 neg=0; end
if nargin<5, xoff=0; end
if nargin<4 phi=0; end
x1=(1:nx)-1;
y1=(1:ny)-1;
[x,y]=meshgrid(x1,y1);
%% the y axis in octave points downward
xp=x*cosd(phi)-y*sind(phi)+xoff;
b=floor(xp/d);
g=0+abs(rem(b,2))==1;
if neg ==0
    g=not(g);
end
endfunction
