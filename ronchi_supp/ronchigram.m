## Copyright (C) 2020 Khosrow Hassani
## University of Tehran, Iran.
## This program is free software; you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 3 of the License, or
## (at your option) any later version.
## 
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*- 
## @deftypefn {} {@var{ga} =} ronchigram (@var{radius}, @var{focal},@var{d}, @var{zoff}, @var{sph}, @var{coma}, @var{astig}, @var{def}, @var{phi}, @var{xoff}, @var{xsize}, @var{ysize}, @var{neg})
##
## Generate Ronchigrams obtained in the Ronchi tests (see "Optical Shop
## testing", by Malacara, ch. 9). RADIUS and FOCAL are the lenns radius
## and focal lenght. D is pitch of the Ronchi ruling lines, and ZOFF is the 
## axial offset of the graing from the paraxial focus (negative toward the lens). SPH, COMA, ASTIG, and DEF are the
## coeficients for spherical, coma, astigmatism, and defocus aberrations (default
## values are zero). PHI is the angle of rotation of the Ronchi grating
## in degrees (default is zero). XSIZE and YSIZE are the X and Y sizes of 
## the returned squared array in pixels (default iz 100 by 100). 
## XOFF is the offset of the center of the Ronchi grating (default is such
## that the central band is a white one). 
## If given and non-zero, NEG will return the negative image (black and white bands switched).
## 
## @seealso{ronchi}
## @end deftypefn

## Author: Khosrow Hassani <hassanikh@ut.ac.ir>
## Created: 2020-10-15

function [ga] = ronchigram (radius,focal,d,zoff,sph,coma,astig,def,phi,xoff,xsize,ysize,neg)
% default values
if nargin<13, neg=0; end
if nargin<12, ysize=100; end 
if nargin<11, xsize=100; end
if nargin<10, xoff=0; end
if nargin<9, phi=0; end
if nargin<8, def=0; end
if nargin<7, astig=0; end
if nargin<6, coma=0; end
if nargin<5, sph=0; end
radius=double(radius);
% define the normalized (to radius) variables:
fn=focal/radius;
zn=zoff/radius;
dn=d/radius;
% image at 2f
rn=2*fn;
% calculate the magnification:
if (zn==0) 
    mag =-rn;
else mag=-rn/zn;
end
% (normalized) exit pupil coordinates:
x1=linspace(-1,1,xsize);
y1=linspace(-1,1,ysize);
[x,y]=meshgrid(x1,y1);
% calculate the aberrated (ray aberration) coordinates:
xp=x-mag*rn*(4*sph*x.*(x.^2+y.^2)+2*coma*x.*y+2*astig*x+2*def*x);
yp=y-mag*rn*(4*sph*y.*(x.^2+y.^2)+coma*x.^2+3*coma*y.^2+6*astig*y+2*def*y);
% define Ronchi ruling pitch in pixels:
d1=round(max(xsize*dn*abs(mag),1));
% set a maximum (normalized) aberrated coordinates big enough
% to accomodate the aberrated coordinates
c=2;
x1=linspace(-c,c,c*xsize);
y1=linspace(-c,c,c*ysize);
[x,y]=meshgrid(x1,y1);
% xoff in pixels
o1=c*xsize*mag*xoff;
% generate the Ronchi ruling on the exit pupil:
g1=ronchi(c*xsize,c*ysize,d1,phi,o1,neg);
% get the aberrated Ronchi ruling (Ronchigram):
rg=interp2(x,y,g1+0.,xp,yp);
% create a circular mask (lens aperture);
x1=linspace(-1,1,xsize);
y1=linspace(-1,1,ysize);
[x,y]=meshgrid(x1,y1);
r=sqrt(x.^2+y.^2);
mask=r<=1;
% mask the Ronchigram
ga=rg.*mask;
endfunction
