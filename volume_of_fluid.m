% =============================================== START =============================================== %

clear; close all; clc;

% [model]
% ^ y
% |
% *-------------*
% |      2      |
% |             |
% |    *---*    |
% |    | 1 |    |
% |    *---*    |
% *-------------*---> x

% [piecewise linear interface calculation]
%   ^ y
%   |
% H *
%   |\ G  C
% D *-*--*
%   |  \2|    nx * x + ny * y = alpha
%   |   \|
%   | 1  * F
%   |    |\
%   *----*-*---> x
%   A    B E

% =========================================== INITIALIZATION ========================================== %

[dt, t] = deal(0.00005, 1);					% time step and domain
[dl, l1, l2] = deal(0.004, 1, 0.4);			% space step and domain
[x_center, y_center] = deal(0.5, 0.3);		% center of phase 1
[left, right, down, up] = deal(x_center - l2/2, x_center + l2/2, y_center - l2/2, y_center + l2/2);	% initial interface

cell_num = l1/dl;
node_num = cell_num + 1;
con = zeros(cell_num);						% concentration field
[conL, conR] = deal(zeros(cell_num));
[conD, conU] = deal(zeros(cell_num));
con(round(down/dl) + 1 : round(up/dl), round(left/dl + 1) : round(right/dl)) = ones(l2/dl);

[x, y] = meshgrid(dl/2 : dl : l1 - dl/2);
[x0, y0] = deal(x - dl/2, y - dl/2);
[x1, y1] = deal(x + dl/2, y + dl/2);
u0 = -2*pi*cos(pi*(x0 - 0.5)).*sin(pi*(y - 0.5));			% velocity field
u1 = -2*pi*cos(pi*(x1 - 0.5)).*sin(pi*(y - 0.5));
v0 =  2*pi*sin(pi*(x - 0.5)).*cos(pi*(y0 - 0.5));
v1 =  2*pi*sin(pi*(x - 0.5)).*cos(pi*(y1 - 0.5));

[node_nx, node_ny] = deal(zeros(node_num));
[nx, ny] = deal(zeros(cell_num));
[Nx, Ny] = deal(zeros(cell_num));
[point1_x, point1_y] = deal(zeros(cell_num));
[point2_x, point2_y] = deal(zeros(cell_num));
alpha = zeros(cell_num);
quad = zeros(cell_num);
many_con = cell(round(t/dt) + 1, 1);
many_con{1} = con;
vol = zeros(round(t/dt) + 1, 1);
vol(1) = sum(con(:))/(cell_num^2);

% ============================================ CALCULATION ============================================ %

for time = 0:dt:t

% -------------------------------------- INTERFACE RECONSTRUCTION ------------------------------------- %

	for node_i = 1:node_num
		for node_j = 1:node_num
			[node_i1, node_j1, node_i2, node_j2] = deal(node_i - 1, node_j - 1, node_i, node_j);
			switch node_i
				case 1
					node_i1 = node_num - 1;
				case node_num
					node_i2 = 1;
			end
			switch node_j
				case 1
					node_j1 = node_num - 1;
				case node_num
					node_j2 = 1;
			end
			node_nx(node_j, node_i) = (con(node_j1, node_i1) + con(node_j2, node_i1) - ...
									   con(node_j2, node_i2) - con(node_j1, node_i2))/2/dl;
			node_ny(node_j, node_i) = (con(node_j1, node_i1) + con(node_j1, node_i2) - ...
									   con(node_j2, node_i2) - con(node_j2, node_i1))/2/dl;
		end
	end

	for i = 1:cell_num
		for j = 1:cell_num
			nx(j, i) = (node_nx(j, i) + node_nx(j + 1, i) + node_nx(j, i + 1) + node_nx(j + 1, i + 1))/4;
			ny(j, i) = (node_ny(j, i) + node_ny(j + 1, i) + node_ny(j, i + 1) + node_ny(j + 1, i + 1))/4;
			[Nx(j, i), Ny(j, i)] = deal(abs(nx(j, i)), abs(ny(j, i)));
		end
	end

	for i = 1:cell_num
		for j = 1:cell_num
			if 	(con(j, i) == 0) || (con(j, i) == 1)
				alpha(j, i) = 0;
			elseif (nx(j, i) == 0) && (ny(j, i) ~= 0)
				alpha(j, i) = con(j, i)*dl*Ny(j, i);
			elseif (ny(j, i) == 0) && (nx(j, i) ~= 0)
				alpha(j, i) = con(j, i)*dl*Nx(j, i);
			elseif (ny(j, i) ~= 0) && (nx(j, i) ~= 0)
				switch sign(Nx(j, i) - Ny(j, i))
					case 0
						if con(j, i) > 0.5
							[DGH, BEF] = deal(1);
						else
							[DGH, BEF] = deal(0);
						end
					case 1
						if con(j, i) <= 0.5*Ny(j, i)/Nx(j, i)
							[DGH, BEF] = deal(0);
						elseif con(j, i) <= 0.5*Nx(j, i)/Ny(j, i)
							[DGH, BEF] = deal(1, 0);
						else
							[DGH, BEF] = deal(1);
						end
					case -1
						if con(j, i) <= 0.5*Nx(j, i)/Ny(j, i)
							[DGH, BEF] = deal(0);
						elseif con(j, i) <= 0.5*Ny(j, i)/Nx(j, i)
							[DGH, BEF] = deal(0, 1);
						else
							[DGH, BEF] = deal(1);
						end
				end

				a = (1 - DGH - BEF)/Nx(j, i)/Ny(j, i);
				b = (DGH/Nx(j, i) + BEF/Ny(j, i))*dl;
				c = -(Ny(j, i)*DGH/2/Nx(j, i) + Nx(j, i)*BEF/2/Ny(j, i) + con(j, i))*(dl^2);
				if a == 0
					alpha(j, i) = -c/b;
				elseif 1 - DGH - BEF == -1
					alpha(j, i) = (Nx(j, i) + Ny(j, i))*dl - sqrt(2*Nx(j, i)*Ny(j, i)*(1 - con(j, i)))*dl;
				else
					alpha(j, i) = sqrt(2*con(j, i)*Nx(j, i)*Ny(j, i))*dl;
				end

				point1_x(j, i) = max(0, (alpha(j, i) - Ny(j, i)*dl)/Nx(j, i));
				point1_y(j, i) = min(alpha(j, i)/Ny(j, i), dl);
				point2_x(j, i) = min(alpha(j, i)/Nx(j, i), dl);
				point2_y(j, i) = max(0, (alpha(j, i) - Nx(j, i)*dl)/Ny(j, i));

				if (nx(j, i) > 0) && (ny(j, i) > 0)
					quad(j, i) = 1;
				elseif (nx(j, i) < 0) && (ny(j, i) > 0)
					quad(j, i) = 2;
				elseif (nx(j, i) < 0) && (ny(j, i) < 0)
					quad(j, i) = 3;
				elseif (nx(j, i) > 0) && (ny(j, i) < 0)
					quad(j, i) = 4;
				end
			end
		end
	end

% -------------------------------------- X-INTERFACE PROPAGATION -------------------------------------- %

	mirror = [2, 3];
	for i = 1:cell_num
		for j = 1:cell_num
			if con(j, i) == 0
				[conL(j, i), conR(j, i)] = deal(0);
			elseif con(j, i) == 1
				conL(j, i) = max(-u0(j, i)*dt/dl, 0);
				conR(j, i) = max(u1(j, i)*dt/dl, 0);
			else
				if (nx(j, i) == 0) && (ny(j, i) ~= 0)
					conL(j, i) = max(-u0(j, i)*dt/dl, 0)*alpha(j, i)/Ny(j, i)/dl;
					conR(j, i) = max(u1(j, i)*dt/dl, 0)*alpha(j, i)/Ny(j, i)/dl;
				elseif (nx(j, i) ~= 0) && (ny(j, i) == 0)
					if nx(j, i) > 0
						ui = (1 - alpha(j, i)/Nx(j, i)/dl)*u0(j, i) + alpha(j, i)/Nx(j, i)/dl*u1(j, i);
						conL(j, i) = max(-u0(j, i)*dt/dl, 0);
						conR(j, i) = max((ui*dt - (dl - alpha(j, i)/Nx(j, i)))/dl, 0);
					else
						ui = (1 - alpha(j, i)/Nx(j, i)/dl)*u1(j, i) + alpha(j, i)/Nx(j, i)/dl*u0(j, i);
						conL(j, i) = max((-ui*dt - (dl - alpha(j, i)/Nx(j, i)))/dl, 0);
						conR(j, i) = max(u1(j, i)*dt/dl, 0);
					end
				else
					if ismember(quad(j, i), mirror)
						[U0, U1] = deal(-u1(j, i), -u0(j, i));
					else
						[U0, U1] = deal(u0(j, i), u1(j, i));
					end
					point1_u = (1 - point1_x(j, i)/dl)*U0 + point1_x(j, i)/dl*U1;
					point2_u = (1 - point2_x(j, i)/dl)*U0 + point2_x(j, i)/dl*U1;
					new_point1_x = point1_x(j, i) + point1_u*dt;
					new_point2_x = point2_x(j, i) + point2_u*dt;
					new_nx = nx(j, i)/(1 + (U1 - U0)*dt/dl);
					new_Nx = abs(new_nx);
					new_alpha = alpha(j, i) + new_nx*U0*dt;
					point3_y = min(new_alpha/Ny(j, i), dl);
					point4_y = max((new_alpha - new_Nx*dl)/Ny(j, i), 0);
					if (new_point1_x < 0) && (new_point2_x < 0)
						conL(j, i) = abs((new_point2_x - U0*dt)*point1_y(j, i)) - ...
									 abs((new_point2_x - new_point1_x)*(point1_y(j, i) - point2_y(j, i)))/2;
						conL(j, i) = conL(j, i)/(dl^2);
						conR(j, i) = 0;
					elseif (new_point1_x > dl) && (new_point2_x > dl)
						conL(j, i) = max(-U0*dt/dl, 0);
						conR(j, i) = abs((new_point2_x - new_point1_x)*(point1_y(j, i) + point2_y(j, i)))/2 + ...
									 min(new_point1_x - dl, abs(new_point1_x - U0*dt))*point1_y(j, i);
						conR(j, i) = conR(j, i)/(dl^2);
					else
						conL(j, i) = max(-U0*dt*point1_y(j, i), 0) - ...
									 max(-abs(point1_y(j, i) - point3_y)*new_point1_x/2, 0);
						conL(j, i) = conL(j, i)/(dl^2);
						conR(j, i) = max((point2_y(j, i) + point4_y)*(new_point2_x - dl)/2, 0);
						conR(j, i) = conR(j, i)/(dl^2);
					end
					if ismember(quad(j, i), mirror)
						[conL(j, i), conR(j, i)] = deal(conR(j, i), conL(j, i));
					end
				end
			end
		end
	end

	for i = 1:cell_num
		[iL, iR] = deal(i - 1, i + 1);
		if i == 1
			iL = cell_num;
		elseif i == cell_num
			iR = 1;
		end
		for j = 1:cell_num
			con(j, i) = con(j, i) - conL(j, i) - conR(j, i) + conL(j, iR) + conR(j, iL);
			if con(j, i) > 0.997
				con(j, i) = 1;
			elseif con(j, i) < 0
				con(j, i) = 0;
			end
		end
	end

% -------------------------------------- INTERFACE RECONSTRUCTION ------------------------------------- %

	for node_i = 1:node_num
		for node_j = 1:node_num
			[node_i1, node_j1, node_i2, node_j2] = deal(node_i - 1, node_j - 1, node_i, node_j);
			switch node_i
				case 1
					node_i1 = node_num - 1;
				case node_num
					node_i2 = 1;
			end
			switch node_j
				case 1
					node_j1 = node_num - 1;
				case node_num
					node_j2 = 1;
			end
			node_nx(node_j, node_i) = (con(node_j1, node_i1) + con(node_j2, node_i1) - ...
									   con(node_j2, node_i2) - con(node_j1, node_i2))/2/dl;
			node_ny(node_j, node_i) = (con(node_j1, node_i1) + con(node_j1, node_i2) - ...
									   con(node_j2, node_i2) - con(node_j2, node_i1))/2/dl;
		end
	end

	for i = 1:cell_num
		for j = 1:cell_num
			nx(j, i) = (node_nx(j, i) + node_nx(j + 1, i) + node_nx(j, i + 1) + node_nx(j + 1, i + 1))/4;
			ny(j, i) = (node_ny(j, i) + node_ny(j + 1, i) + node_ny(j, i + 1) + node_ny(j + 1, i + 1))/4;
			[Nx(j, i), Ny(j, i)] = deal(abs(nx(j, i)), abs(ny(j, i)));
		end
	end

	for i = 1:cell_num
		for j = 1:cell_num
			if 	(con(j, i) == 0) || (con(j, i) == 1)
				alpha(j, i) = 0;
			elseif (nx(j, i) == 0) && (ny(j, i) ~= 0)
				alpha(j, i) = con(j, i)*dl*Ny(j, i);
			elseif (ny(j, i) == 0) && (nx(j, i) ~= 0)
				alpha(j, i) = con(j, i)*dl*Nx(j, i);
			elseif (ny(j, i) ~= 0) && (nx(j, i) ~= 0)
				switch sign(Nx(j, i) - Ny(j, i))
					case 0
						if con(j, i) > 0.5
							[DGH, BEF] = deal(1);
						else
							[DGH, BEF] = deal(0);
						end
					case 1
						if con(j, i) <= 0.5*Ny(j, i)/Nx(j, i)
							[DGH, BEF] = deal(0);
						elseif con(j, i) <= 0.5*Nx(j, i)/Ny(j, i)
							[DGH, BEF] = deal(1, 0);
						else
							[DGH, BEF] = deal(1);
						end
					case -1
						if con(j, i) <= 0.5*Nx(j, i)/Ny(j, i)
							[DGH, BEF] = deal(0);
						elseif con(j, i) <= 0.5*Ny(j, i)/Nx(j, i)
							[DGH, BEF] = deal(0, 1);
						else
							[DGH, BEF] = deal(1);
						end
				end

				a = (1 - DGH - BEF)/Nx(j, i)/Ny(j, i);
				b = (DGH/Nx(j, i) + BEF/Ny(j, i))*dl;
				c = -(Ny(j, i)*DGH/2/Nx(j, i) + Nx(j, i)*BEF/2/Ny(j, i) + con(j, i))*(dl^2);
				if a == 0
					alpha(j, i) = -c/b;
				elseif 1 - DGH - BEF == -1
					alpha(j, i) = (Nx(j, i) + Ny(j, i))*dl - sqrt(2*Nx(j, i)*Ny(j, i)*(1 - con(j, i)))*dl;
				else
					alpha(j, i) = sqrt(2*con(j, i)*Nx(j, i)*Ny(j, i))*dl;
				end

				point1_x(j, i) = max(0, (alpha(j, i) - Ny(j, i)*dl)/Nx(j, i));
				point1_y(j, i) = min(alpha(j, i)/Ny(j, i), dl);
				point2_x(j, i) = min(alpha(j, i)/Nx(j, i), dl);
				point2_y(j, i) = max(0, (alpha(j, i) - Nx(j, i)*dl)/Ny(j, i));

				if (nx(j, i) > 0) && (ny(j, i) > 0)
					quad(j, i) = 1;
				elseif (nx(j, i) < 0) && (ny(j, i) > 0)
					quad(j, i) = 2;
				elseif (nx(j, i) < 0) && (ny(j, i) < 0)
					quad(j, i) = 3;
				elseif (nx(j, i) > 0) && (ny(j, i) < 0)
					quad(j, i) = 4;
				end
			end
		end
	end

% -------------------------------------- Y-INTERFACE PROPAGATION -------------------------------------- %

	mirror = [3, 4];
	for i = 1:cell_num
		for j = 1:cell_num
			if con(j, i) == 0
				[conD(j, i), conU(j, i)] = deal(0);
			elseif con(j, i) == 1
				conD(j, i) = max(-v0(j, i)*dt/dl, 0);
				conU(j, i) = max(v1(j, i)*dt/dl, 0);
			else
				if (nx(j, i) == 0) && (ny(j, i) ~= 0)
					if ny(j, i) > 0
						vi = (1 - alpha(j, i)/Ny(j, i)/dl)*v0(j, i) + alpha(j, i)/Ny(j, i)/dl*v1(j, i);
						conD(j, i) = max(-v0(j, i)*dt/dl, 0);
						conU(j, i) = max((vi*dt - (dl - alpha(j, i)/Ny(j, i)))/dl, 0);
					else
						vi = (1 - alpha(j, i)/Ny(j, i)/dl)*v1(j, i) + alpha(j, i)/Ny(j, i)/dl*v0(j, i);
						conD(j, i) = max((-vi*dt - (dl - alpha(j, i)/Ny(j, i)))/dl, 0);
						conU(j, i) = max(v1(j, i)*dt/dl, 0);
					end
				elseif (nx(j, i) ~= 0) && (ny(j, i) == 0)
					conD(j, i) = max(-v0(j, i)*dt/dl, 0)*alpha(j, i)/Nx(j, i)/dl;
					conU(j, i) = max(v1(j, i)*dt/dl, 0)*alpha(j, i)/Nx(j, i)/dl;
				else
					if ismember(quad(j, i), mirror)
						[V0, V1] = deal(-v1(j, i), -v0(j, i));
					else
						[V0, V1] = deal(v0(j, i), v1(j, i));
					end
					point1_v = (1 - point1_y(j, i)/dl)*V0 + point1_y(j, i)/dl*V1;
					point2_v = (1 - point2_y(j, i)/dl)*V0 + point2_y(j, i)/dl*V1;
					new_point1_y = point1_y(j, i) + point1_v*dt;
					new_point2_y = point2_y(j, i) + point2_v*dt;
					new_ny = ny(j, i)/(1 + (V1 - V0)*dt/dl);
					new_Ny = abs(new_ny);
					new_alpha = alpha(j, i) + new_ny*V0*dt;
					point3_x = max((new_alpha - new_Ny*dl)/Nx(j, i), 0);
					point4_x = min(new_alpha/Nx(j, i), dl);
					if (new_point1_y < 0) && (new_point2_y < 0)
						conD(j, i) = abs((new_point1_y - V0*dt)*point2_x(j, i)) - ...
									 abs((new_point1_y - new_point2_y)*(point2_x(j, i) - point1_x(j, i)))/2;
						conD(j, i) = conD(j, i)/(dl^2);
						conU(j, i) = 0;
					elseif (new_point1_y > dl) && (new_point2_y > dl)
						conD(j, i) = max(-V0*dt/dl, 0);
						conU(j, i) = abs((point1_x(j, i) + point2_x(j, i))*(new_point1_y - new_point2_y))/2 + ...
									 min(new_point2_y - dl, abs(new_point2_y - V0*dt))*point2_x(j, i);
						conU(j, i) = conU(j, i)/(dl^2);
					else
						conD(j, i) = max(-V0*dt*point2_x(j, i), 0) - ...
									 max(-abs(point2_x(j, i) - point4_x)*new_point2_y/2, 0);
						conD(j, i) = conD(j, i)/(dl^2);
						conU(j, i) = max((point1_x(j, i) + point3_x)*(new_point1_y - dl)/2, 0);
						conU(j, i) = conU(j, i)/(dl^2);
					end
					if ismember(quad(j, i), mirror)
						[conD(j, i), conU(j, i)] = deal(conU(j, i), conD(j, i));
					end
				end
			end
		end
	end

	for i = 1:cell_num
		for j = 1:cell_num
			[jD, jU] = deal(j - 1, j + 1);
			if j == 1
				jD = cell_num;
			elseif j == cell_num
				jU = 1;
			end
			con(j, i) = con(j, i) - conD(j, i) - conU(j, i) + conD(jU, i) + conU(jD, i);
            if con(j, i) > 0.997
				con(j, i) = 1;
			elseif con(j, i) < 0
				con(j, i) = 0;
            end
		end
	end

	if time > 0
		many_con{1 + round(time/dt)} = con;
	end
	vol(1 + round(time/dt)) = sum(con(:))/(cell_num^2);
	disp(['step = ', num2str(round(time/dt) + 1), '; ', ...
		  'time = ', num2str(time, '%6.5f'), '; ', ...
		  'volume = ', num2str(vol(1 + round(time/dt)), '%7.6f')]);
end

disp('ALL OVER!');
% ================================================ END ================================================ %