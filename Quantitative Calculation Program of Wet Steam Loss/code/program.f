!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	program wetness loss
	use WASPCN
	implicit double precision (a-h,o-z)
	parameter (zhuansu=1500*2*3.1415926/60.0)
!	zhuansu:ת��,���ٶ�
	parameter (NStage=10,NSpan=11,NStream=11)
!	NStage:����,NSpan:չ��ڵ���,NStream:��Ҷդ����ڵ���
	parameter (NNode=2*NStage*NSpan*NStream)
	parameter (MNode=2*NStage*NStream)
	double precision MassFlow
c	parameter (MassFlow=104.03412)
	parameter (MassFlow=271.8685)
	double precision NIndex(0:NNode-1)
!	NNode:�ܽڵ���Ŀ
	parameter (MS=100,MN=800)
!	���������������������Ƶ����
	double precision x(0:NNode-1),y(0:NNode-1),z(0:NNode-1)
	double precision u_a(0:NNode-1),u_c(0:NNode-1),u_r(0:NNode-1)
	double precision rho_v(0:NNode-1),rho_l(0:NNode-1)
	double precision p_d(0:NNode-1),d_v(0:NNode-1),d_n(0:NNode-1)
	double precision v_f(0:NNode-1),s_s(0:NNode-1),s_t(0:NNode-1)
	double precision t_k_e(0:NNode-1)
!	�����ʼ�����ݶ���
	double precision xx(0:NNode-1),yy(0:NNode-1),zz(0:NNode-1)
	double precision uu_axial(0:NNode-1),uu_circumferential(0:NNode-1)
	double precision uu_radial(0:NNode-1)
	double precision rho_vv(0:NNode-1),rho_ll(0:NNode-1)
	double precision p_dd(0:NNode-1),d_vv(0:NNode-1),d_nn(0:NNode-1)
	double precision v_ff(0:NNode-1),s_ss(0:NNode-1),s_tt(0:NNode-1)
	double precision t_kk_ee(0:NNode-1)
!	����˳����ʱ���ݶ���
	common/com01/cylindrical_r(0:2*NStage*NStream-1,0:NSpan-1)
	common/com02/cylindrical_z(0:2*NStage*NStream-1,0:NSpan-1)
	common/com03/omiga(0:2*NStage*NStream-1,0:NSpan-1)
	common/com04/velocity_axial(0:2*NStage*NStream-1,0:NSpan-1)
	common/com05/velocity_radial(0:2*NStage*NStream-1,0:NSpan-1)
	common/com06/velocity_circumferential(0:2*NStage*NStream-1,
     &	0:NSpan-1)
	common/com07/density_vapor(0:2*NStage*NStream-1,0:NSpan-1)
	common/com08/density_liquid(0:2*NStage*NStream-1,0:NSpan-1)
	common/com09/particle_diameter(0:2*NStage*NStream-1,0:NSpan-1)
	common/com10/dynamic_viscosity(0:2*NStage*NStream-1,0:NSpan-1)
	common/com11/droplet_number(0:2*NStage*NStream-1,0:NSpan-1)
	common/com12/streamwise_location(0:2*NStage*NStream-1)
	common/com13/static_pressure(0:2*NStage*NStream-1)
	common/com14/static_temperature(0:2*NStage*NStream-1)
	common/com15/d_mass_fraction(0:2*NStage*NStream-1)
	common/com16/axial_distance(0:2*NStage*NStream-1)
	common/com17/a_mass_flowrate(0:2*NStage*NStream-1)
	common/com18/volume_fraction(0:2*NStage*NStream-1,0:NSpan-1)
	common/com19/sound_speed(0:2*NStage*NStream-1,0:NSpan-1)
	common/com20/surface_tension(0:2*NStage*NStream-1,0:NSpan-1)
	common/com21/turbulence_kinetic_energy(0:2*NStage*NStream-1,
     &	0:NSpan-1)
!	����˳���ṹ�����ݶ���
	common/com30/v_d_a(0:2*NStage*NStream-1,0:NSpan-1)
	common/com31/v_d_r(0:2*NStage*NStream-1,0:NSpan-1)
	common/com32/v_d_c(0:2*NStage*NStream-1,0:NSpan-1)
!	һ��Һ��������ϵ���������ٶ�
	integer bd(0:NSpan-1)
!	���������Ὺʼ������λ��
	common/com33/v_ds_a(0:2*NStage*NStream-1,0:NSpan-1)
	common/com34/v_ds_r(0:2*NStage*NStream-1,0:NSpan-1)
	common/com35/v_ds_c(0:2*NStage*NStream-1,0:NSpan-1)
!	һ��Һ��������ϵ���������ٶ�
	double precision k1,k2,k3,k4
	double precision v_l_r(0:1000)
	double precision v_l_c(0:1000)
	double precision v_l_a(0:1000000)
	double precision p_l_a(0:1000000)
!	�����򵥶����ٶȺ�λ����ʱ����a:����,r:����,c:����
	double precision qx_weizhi(0:NSpan*NStream-1)
	double precision qx_sudu(0:NSpan*NStream-1)
	double precision qx_midu(0:NSpan*NStream-1)
	double precision qx_nianxing(0:NSpan*NStream-1)
	double precision qx_jiaosudu(0:NSpan*NStream-1)
	double precision kl_midu(0:NSpan*NStream-1)
	double precision kl_zhijing(0:NSpan*NStream-1)
	double precision kl_c_sudu(0:NSpan*NStream-1)
	double precision kl_surfacetension(0:NSpan*NStream-1)
	double precision kl_shengsu(0:NSpan*NStream-1)
	double precision kl_tijifenshu(0:NSpan*NStream-1)
	double precision qxsudu
	double precision qxmidu
	double precision qxnianxing
	double precision qxweizhi
	double precision qxjiaosudu
	double precision klmidu
	double precision klzhijing
	double precision klcsudu
!	����Һ����������򵥶�����ʱ�ٶȺ�λ��
	common/com40/Number_s(NStage)
	common/com41/Number_r(NStage)
!	������Ҷ�Ͷ�ҶҶƬ����
	double precision chushisudu(NStage,0:NSpan)
	double precision chushiweizhi(NStage,0:NSpan)
	double precision chushiqxmidu
	double precision chushiqxsudu
	double precision klsurfacetension
	double precision kldmidu(NStage,0:NSpan-1)
	double precision kldzhijing(NStage,0:NSpan-1)
	double precision klsezhijing
	double precision klsemidu
	double precision zhidongjiao(NStage,0:NSpan)
	double precision temp
	double precision tmp(2*NStage)
	double precision xs(0:NNode-1),ys(0:NNode-1),zs(0:NNode-1)
	double precision xr(0:NNode-1),yr(0:NNode-1),zr(0:NNode-1)
!	��Ҷ�Ͷ�ҶҶ������
	character*50 line
	common/com50/NSection_s(NStage)
	common/com51/NSection_r(NStage)
!	��Ҷ�����������
	common/com52/NSNumber_s(NStage,0:MS)
	common/com53/NPNumber_s(NStage,0:MS)
	common/com54/NSNumber_r(NStage,0:MS)
	common/com55/NPNumber_r(NStage,0:MS)
!	��Ҷ��������Ƶ��������Ҷѹ������Ƶ����
	common/com57/NSNumber_rn(NStage,0:MS)
!	������Ҷ�������޳��ǵ������к�Ŀ��Ƶ����
	integer NSNumber_rns(NStage)
!	������Ҷ��������С���Ƶ����
	character abc *20, def *20
	character filename *20 (2*NStage)
	integer aaa
!	�ַ���:�ļ���
	common/com56/idirection
!	��Ҷ��ת�������жϱ���
	double precision xss(NStage,0:MS,0:MN)
	double precision yss(NStage,0:MS,0:MN)
	double precision zss(NStage,0:MS,0:MN)
	double precision xps(NStage,0:MS,0:MN)
	double precision yps(NStage,0:MS,0:MN)
	double precision zps(NStage,0:MS,0:MN)
!	��Ҷ�������ѹ����ֱ������������������Ƶ�����
	double precision xsr(NStage,0:MS,0:MN)
	double precision ysr(NStage,0:MS,0:MN)
	double precision zsr(NStage,0:MS,0:MN)
	double precision xpr(NStage,0:MS,0:MN)
	double precision ypr(NStage,0:MS,0:MN)
	double precision zpr(NStage,0:MS,0:MN)
!	��Ҷ�������ѹ����ֱ������������������Ƶ�����
	common/com60/rss(NStage,0:MS,0:MN)
	common/com61/tss(NStage,0:MS,0:MN)
	common/com62/ass(NStage,0:MS,0:MN)
	common/com63/rps(NStage,0:MS,0:MN)
	common/com64/tps(NStage,0:MS,0:MN)
	common/com65/aps(NStage,0:MS,0:MN)
!	��Ҷ�������ѹ����Բ��������Ƶ�����
	common/com66/rsr(NStage,0:MS,0:MN)
	common/com67/tsr(NStage,0:MS,0:MN)
	common/com68/asr(NStage,0:MS,0:MN)
	common/com69/rpr(NStage,0:MS,0:MN)
	common/com70/tpr(NStage,0:MS,0:MN)
	common/com71/apr(NStage,0:MS,0:MN)
!	��Ҷ�������ѹ����Բ��������Ƶ�����
	integer istart,iend,ns
	common/com72/rsrn(NStage,0:MS,0:MN)
	common/com73/tsrn(NStage,0:MS,0:MN)
	common/com74/asrn(NStage,0:MS,0:MN)
!	��Ҷ������Բ�������޳�ǰ��Ե����ǵ������Ƶ�����
	double precision rr(NStage,0:NSpan-1,0:MN)
	double precision rz(NStage,0:NSpan-1,0:MN)
	double precision rc(NStage,0:NSpan-1,0:MN)
!	��Ҷ�������ض�Ҷ�߾���,�������������
	double precision slr(NStage,0:NSpan-1),slz(NStage,0:NSpan-1)
	double precision str(NStage,0:NSpan-1),stz(NStage,0:NSpan-1)
	double precision rlr(NStage,0:NSpan-1),rlz(NStage,0:NSpan-1)
	double precision rtr(NStage,0:NSpan-1),rtz(NStage,0:NSpan-1)
!	��Ҷ�Ͷ�Ҷǰ��Ե�������������
	common/com80/turbine_power(NStage)
!	��������(�ֶ��������)
	double precision f_v(0:2*NStage*NStream-1,0:NSpan-1)
!	�����ٶ�
	double precision r,c,za
!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	������������
!	���������ٶ�
	open(1,file='Velocity Axial.csv')
	do i=1,6
		read(1,10) line
	enddo
	do i=0,NNode-1
		read(1,*) NIndex(i),x(i),y(i),z(i),u_a(i)
	enddo
	close(1)
!	���������ٶ�
	open(2,file='Velocity Circumferential.csv')
	do i=1,6
		read(2,10) line
	enddo
	do i=0,NNode-1
		read(2,*) NIndex(i),x(i),y(i),z(i),u_c(i)
	enddo
	close(2)
!	���뾶���ٶ�
	open(3,file='Velocity Radial.csv')
	do i=1,6
		read(3,10) line
	enddo
	do i=0,NNode-1
		read(3,*) NIndex(i),x(i),y(i),z(i),u_r(i)
	enddo
	close(3)
!	���������ܶ�
	open(4,file='Vapor Density.csv')
	do i=1,6
		read(4,10) line
	enddo
	do i=0,NNode-1
		read(4,*) NIndex(i),x(i),y(i),z(i),rho_v(i)
	enddo
	close(4)
!	����Һ���ܶ�
	open(5,file='Liquid Density.csv')
	do i=1,6
		read(5,10) line
	enddo
	do i=0,NNode-1
		read(5,*) NIndex(i),x(i),y(i),z(i),rho_l(i)
	enddo
	close(5)
!	����Һ��ֱ��
	open(6,file='Particle Diameter.csv')
	do i=1,6
		read(6,10) line
	enddo
	do i=0,NNode-1
		read(6,*) NIndex(i),x(i),y(i),z(i),p_d(i)
	enddo
	close(6)
!	�������ද��ճ��
	open(7,file='Dynamic Viscosity.csv')
	do i=1,6
		read(7,10) line
	enddo
	do i=0,NNode-1
		read(7,*) NIndex(i),x(i),y(i),z(i),d_v(i)
	enddo
	close(7)
!	����Һ����Ŀ
	open(8,file='Droplet Number.csv')
	do i=1,6
		read(8,10) line
	enddo
	do i=0,NNode-1
		read(8,*) NIndex(i),x(i),y(i),z(i),d_n(i)
	enddo
	close(8)
!	��������ѹ
	open(9,file='Static Pressure.csv')
	do i=1,5
		read(9,10) line
	enddo
	do i=0,MNode-1
		read(9,*) streamwise_location(i),static_pressure(i)
	enddo
	close(9)
!	����������
	open(10,file='Static Temperature.csv')
	do i=1,5
		read(10,10) line
	enddo
	do i=0,MNode-1
		read(10,*) streamwise_location(i),static_temperature(i)
	enddo
	close(10)
!	��������ʪ��
	open(11,file='Mass Fraction.csv')
	do i=1,5
		read(11,10) line
	enddo
	do i=0,MNode-1
		read(11,*) streamwise_location(i),d_mass_fraction(i)
	enddo
	close(11)
!	������������
	open(12,file='Axial Distance.csv')
	do i=1,5
		read(12,10) line
	enddo
	do i=0,MNode-1
		read(12,*) streamwise_location(i),axial_distance(i)
	enddo
	close(12)
!	����Һ���������
	open(13,file='Volume Fraction.csv')
	do i=1,6
		read(13,10) line
	enddo
	do i=0,NNode-1
		read(13,*) NIndex(i),x(i),y(i),z(i),v_f(i)
	enddo
	close(13)
!	����Һ������
	open(15,file='Sound Speed.csv')
	do i=1,6
		read(15,10) line
	enddo
	do i=0,NNode-1
		read(15,*) NIndex(i),x(i),y(i),z(i),s_s(i)
	enddo
	close(15)
!	�����������
	open(16,file='Surface Tension.csv')
	do i=1,6
		read(16,10) line
	enddo
	do i=0,NNode-1
		read(16,*) NIndex(i),x(i),y(i),z(i),s_t(i)
	enddo
	close(16)
!	�����Ķ���
	open(17,file='Turbulence Kinetic Energy.csv')
	do i=1,6
		read(17,10) line
	enddo
	do i=0,NNode-1
		read(17,*) NIndex(i),x(i),y(i),z(i),t_k_e(i)
	enddo
	close(17)
!	��ȡҶ�ͼ�������
!!!!!!
c	character abc *20, def *20
c	character filename *20 (2*NStage)
c	integer aaa
!	�ַ���:�ļ���
!!!!!!
c	abc='s1.geomTurbo'
	abc='S1.geomTurbo'
	def=abc
	aaa=ichar(abc(2:2))-1
	do k=1,NStage
		i=2*k-1
		j=2*k
c		def(1:1)='d'
		def(1:1)='R'
		abc(2:2)=char(aaa+k)
		def(2:2)=char(aaa+k)
		filename(i)=abc
		filename(j)=def
c	print *,i,filename(i),j,filename(j)
	enddo
!	���ļ����ַ�����ֵ
	filename(19)='S10.geomTurbo'
	filename(20)='R10.geomTurbo'
c	print *,filename
!!!!!!�ļ��ַ���������!
	do k=1,2*NStage
	open(20+k,file=filename(k))
	if (mod(k,2)==1) then
!	����Ҷ����
	do i=1,6
		read(20+k,10)line
	enddo
	read(20+k,*),line,Number_s(k/2+1)
	do i=8,13
		read(20+k,10)line
	enddo
	read(20+k,*)
	read(20+k,*)
	read(20+k,*)NSection_s(k/2+1)
	do i=0,NSection_s(k/2+1)-1
		read(20+k,*)
		read(20+k,*)
		read(20+k,*)NSNumber_s(k/2+1,i)
			do j=0,NSNumber_s(k/2+1,i)-1
			read(20+k,*)xss(k/2+1,i,j),yss(k/2+1,i,j),zss(k/2+1,i,j)
			enddo
	enddo
	read(20+k,*)
	read(20+k,*)
	read(20+k,*)
	do i=0,NSection_s(k/2+1)-1
		read(20+k,*)
		read(20+k,*)
		read(20+k,*)NPNumber_s(k/2+1,i)
			do j=0,NPNumber_s(k/2+1,i)-1
			read(20+k,*)xps(k/2+1,i,j),yps(k/2+1,i,j),zps(k/2+1,i,j)
			enddo
	enddo
	elseif (mod(k,2)==0) then
!	����Ҷ����
	do i=1,6
		read(20+k,10)line
	enddo
	read(20+k,*),line,Number_r(k/2)
	do i=8,13
		read(20+k,10)line
	enddo
	read(20+k,*)
	read(20+k,*)
	read(20+k,*)NSection_r(k/2)
	do i=0,NSection_r(k/2)-1
		read(20+k,*)
		read(20+k,*)
		read(20+k,*)NSNumber_r(k/2,i)
			do j=0,NSNumber_r(k/2,i)-1
			read(20+k,*)xsr(k/2,i,j),ysr(k/2,i,j),zsr(k/2,i,j)
			enddo
	enddo
	read(20+k,*)
	read(20+k,*)
	read(20+k,*)
	do i=0,NSection_r(k/2)-1
		read(20+k,*)
		read(20+k,*)
		read(20+k,*)NPNumber_r(k/2,i)
			do j=0,NPNumber_r(k/2,i)-1
			read(20+k,*)xpr(k/2,i,j),ypr(k/2,i,j),zpr(k/2,i,j)
			enddo
	enddo
	endif
	close(20+k)	
	enddo
10	format(A)
!	�������ݶ������
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	�����м���ֱ�����꼸�β���ת��Ϊ�����꣬����������
!	����һ����Ҷ�Ͷ�Ҷֱ������ת��Ϊ������(��λ����Ϊm)
	do k=1,NStage
	do i=0,NSection_s(k)-1
		do j=0,NSNumber_s(k,i)-1
			rss(k,i,j)=sqrt(xss(k,i,j)**2+yss(k,i,j)**2)/1000.0
	tss(k,i,j)=rss(k,i,j)*acos(xss(k,i,j)/1000.0/rss(k,i,j))
			ass(k,i,j)=zss(k,i,j)/1000.0
		enddo
		do j=0,NPNumber_s(k,i)-1
			rps(k,i,j)=sqrt(xps(k,i,j)**2+yps(k,i,j)**2)/1000.0
	tps(k,i,j)=rps(k,i,j)*acos(xps(k,i,j)/1000.0/rps(k,i,j))
			aps(k,i,j)=zps(k,i,j)/1000.0
		enddo
	enddo
	do i=0,NSection_r(k)-1
		do j=0,NSNumber_r(k,i)-1
			rsr(k,i,j)=sqrt(xsr(k,i,j)**2+ysr(k,i,j)**2)/1000.0
	tsr(k,i,j)=rsr(k,i,j)*acos(xsr(k,i,j)/1000.0/rsr(k,i,j))
			asr(k,i,j)=zsr(k,i,j)/1000.0
		enddo
		do j=0,NPNumber_r(k,i)-1
			rpr(k,i,j)=sqrt(xpr(k,i,j)**2+ypr(k,i,j)**2)/1000.0
	tpr(k,i,j)=rpr(k,i,j)*acos(xpr(k,i,j)/1000.0/rpr(k,i,j))
			apr(k,i,j)=zpr(k,i,j)/1000.0
		enddo
	enddo
	enddo
!!!!!!
!	�Զ�Ҷ���������ݽ���ɾѡ���޳����˵�ǵ�������
	do k=1,NStage
		do i=0,NSection_r(k)-1
			do j=0,NSNumber_r(k,i)-2
				if (asr(k,i,j)<asr(k,i,j+1)) then
				iend=j+1
!	���㵥�������������							
				endif 				
			enddo
			do j=NSNumber_r(k,i)-2,0,-1
				if (asr(k,i,j)<asr(k,i,j+1)) then
				istart=j		
!	���㵥��������С����							
				endif 				
			enddo
			NSNumber_rn(k,i)=iend-istart+1
			do j=0,NSNumber_rn(k,i)-1
				rsrn(k,i,j)=rsr(k,i,istart+j)
				tsrn(k,i,j)=tsr(k,i,istart+j)
				asrn(k,i,j)=asr(k,i,istart+j)
!	���µĶ�Ҷ���������ݸ�ֵ
			enddo
		enddo
	enddo
!!!!!!
!	�ж϶�Ҷ��ת������
	if (tss(1,0,NSNumber_s(1,0)-1)>tss(1,0,0)) then
		idirection=1
!	��ת����������ַ���Ϊ������
	else
		idirection=-1
!	��ת����������ַ���Ϊ������
	endif
c	do i=1,2
c		do j=2*(i-1)*NStream,2*i*NStream-1
c			a_mass_flowrate(j)=MassFlow
c		enddo
c	enddo
c	do i=3,3
c		do j=2*(i-1)*NStream,2*i*NStream-1
c			a_mass_flowrate(j)=MassFlow-3.705
c		enddo
c	enddo
c	do i=4,5
c		do j=2*(i-1)*NStream,2*i*NStream-1
c			a_mass_flowrate(j)=MassFlow-3.705-7.365
c		enddo
c	enddo
c	do i=6,7
c		do j=2*(i-1)*NStream,2*i*NStream-1
c			a_mass_flowrate(j)=MassFlow-3.705-7.365-7.206
c		enddo
c	enddo
!	����������ֵ
	do i=1,3
		do j=2*(i-1)*NStream,2*i*NStream-1
			a_mass_flowrate(j)=MassFlow
		enddo
	enddo
	do i=4,7
		do j=2*(i-1)*NStream,2*i*NStream-1
			a_mass_flowrate(j)=MassFlow-25.5325
		enddo
	enddo
	do i=8,9
		do j=2*(i-1)*NStream,2*i*NStream-1
			a_mass_flowrate(j)=MassFlow-25.5325-18.8270
		enddo
	enddo
	do i=10,10
		do j=2*(i-1)*NStream,2*i*NStream-1
			a_mass_flowrate(j)=MassFlow-25.5325-18.8270-10.7645
		enddo
	enddo
!	����������ֵ
!!!!!!
c	turbine_power(1)=11232168.17
c	turbine_power(2)=11396732.83
c	turbine_power(3)=11309242.62
c	turbine_power(4)=10356854.73
c	turbine_power(5)=10881550.28
c	turbine_power(6)=11599569.33
c	turbine_power(7)=13159951.76
!	�������ʴ�С
	turbine_power(1)=15520800
	turbine_power(2)=15028000
	turbine_power(3)=16116800
	turbine_power(4)=14640200
	turbine_power(5)=15703000
	turbine_power(6)=17602400
	turbine_power(7)=21001500
	turbine_power(8)=20203900
	turbine_power(9)=23321800
	turbine_power(10)=26420700
!	�������ʴ�С
!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	������Ҷ�Ͷ�Ҷ�����������������˳��
	call adjust_data_order(NStage,NSpan,NStream,NNode,z,x,xx)
	call adjust_data_order(NStage,NSpan,NStream,NNode,z,y,yy)
	call adjust_data_order(NStage,NSpan,NStream,NNode,z,z,zz)
	call adjust_data_order(NStage,NSpan,NStream,NNode,z,
     &	u_a,uu_axial)
	call adjust_data_order(NStage,NSpan,NStream,NNode,z,
     &	u_c,uu_circumferential)
	call adjust_data_order(NStage,NSpan,NStream,NNode,z,
     &	u_r,uu_radial)
	call adjust_data_order(NStage,NSpan,NStream,NNode,z,
     &	rho_v,rho_vv)
	call adjust_data_order(NStage,NSpan,NStream,NNode,z,
     &	rho_l,rho_ll)
	call adjust_data_order(NStage,NSpan,NStream,NNode,z,
     &	p_d,p_dd)
	call adjust_data_order(NStage,NSpan,NStream,NNode,z,
     &	d_v,d_vv)
	call adjust_data_order(NStage,NSpan,NStream,NNode,z,
     &	d_n,d_nn)
	call adjust_data_order(NStage,NSpan,NStream,NNode,z,
     &	v_f,v_ff)
	call adjust_data_order(NStage,NSpan,NStream,NNode,z,
     &	s_s,s_ss)
	call adjust_data_order(NStage,NSpan,NStream,NNode,z,
     &	s_t,s_tt)
	call adjust_data_order(NStage,NSpan,NStream,NNode,z,
     &	t_k_e,t_kk_ee)
!	��ֱ������ת��ΪԲ������
	do i=0,2*NStage*NStream-1
		do j=0,NSpan-1
			cylindrical_r(i,j)=
     &		sqrt((xx(NSpan*i+j))**2+(yy(NSpan*i+j))**2)
			cylindrical_z(i,j)=zz(NSpan*i+j)
		enddo
	enddo
!	�������ݸ�ʽת��
	call adjust_coordinate(NStage,NSpan,NStream,NNode,
     &	uu_axial,velocity_axial)
	call adjust_coordinate(NStage,NSpan,NStream,NNode,
     &	uu_circumferential,velocity_circumferential)
	call adjust_coordinate(NStage,NSpan,NStream,NNode,
     &	uu_radial,velocity_radial)
	call adjust_coordinate(NStage,NSpan,NStream,NNode,
     &	rho_vv,density_vapor)
	call adjust_coordinate(NStage,NSpan,NStream,NNode,
     &	rho_ll,density_liquid)
	call adjust_coordinate(NStage,NSpan,NStream,NNode,
     &	p_dd,particle_diameter)
	call adjust_coordinate(NStage,NSpan,NStream,NNode,
     &	d_vv,dynamic_viscosity)
	call adjust_coordinate(NStage,NSpan,NStream,NNode,
     &	d_nn,droplet_number)
	call adjust_coordinate(NStage,NSpan,NStream,NNode,
     &	v_ff,volume_fraction)
	call adjust_coordinate(NStage,NSpan,NStream,NNode,
     &	s_ss,sound_speed)
	call adjust_coordinate(NStage,NSpan,NStream,NNode,
     &	s_tt,surface_tension)
	call adjust_coordinate(NStage,NSpan,NStream,NNode,
     &	t_kk_ee,turbulence_kinetic_energy)
!	���ٶȶ���
	do j=0,NSpan-1
		do k=0,NStage-1
			do i=0,NStream-1
				omiga(2*k*NStream+i,j)=0.0
				omiga((2*k+1)*NStream+i,j)=zhuansu
			enddo
		enddo
	enddo
!	�����������ٶ�ת��Ϊ�����ٶ�
	do i=0,2*NStage*NStream-1
		do j=0,NSpan-1
			velocity_circumferential(i,j)=
     &		velocity_circumferential(i,j)+
     &		cylindrical_r(i,j)*omiga(i,j)
		enddo
	enddo
!	�ڼ�������߽紦����������ͱ���ֵ���д������������������ֵ��ͬ
	call DATA_ADJUST_COORDINATE(NStage,NSpan,NStream,cylindrical_z)
	call DATA_ADJUST_VALUE(NStage,NSpan,NStream,density_vapor)
	call DATA_ADJUST_VALUE(NStage,NSpan,NStream,particle_diameter)
	call DATA_ADJUST_VALUE(NStage,NSpan,NStream,dynamic_viscosity)
	call DATA_ADJUST_VALUE(NStage,NSpan,NStream,volume_fraction)
	call DATA_ADJUST_VALUE(NStage,NSpan,NStream,cylindrical_r)

	call DATA_ADJUST_VALUE(NStage,NSpan,NStream,velocity_radial)
	call DATA_ADJUST_VALUE(NStage,NSpan,NStream,
     &	velocity_circumferential)
!	call DATA_ADJUST_VALUE(NStage,NSpan,NStream,velocity_axial)
!	�������������Ὺʼ������λ��
	call condensate_start(bd)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	Һ��������ʼ�ٶ�
	do j=0,NSpan-1
		do i=0,2*NStage*NStream-1
			v_d_a(i,j)=velocity_axial(i,j)
			v_d_r(i,j)=velocity_radial(i,j)
			v_d_c(i,j)=velocity_circumferential(i,j)
		enddo
	enddo
!	����Һ��������ʼ�ٶ�
	do j=0,NSpan-1
		do i=0,2*NStage*NStream-1
			v_ds_a(i,j)=0.0
			v_ds_r(i,j)=0.0
			v_ds_c(i,j)=0.0
		enddo
	enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	����Runge Kutta������һ��Һ���ٶȳ�
	print *,'��ʼ����һ��Һ���ٶ�......'
!	���������ٶ�
	call fluctuating_velocity(f_v)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	���㾶��һ��Һ���ٶȳ�
	do 300,j=0,NSpan-1
	do 250,i=bd(j)+1,2*NStage*NStream-1
	v_l_r(0)=velocity_radial(i,j)
	v_l_r(1)=velocity_radial(i,j)+0.01
!	�����ʼ�ٶ�
	kk=1
	do 270 while (abs(v_l_r(kk)-v_l_r(kk-1))>1D-6)
	k1=delta_r(density_vapor(i,j),(velocity_radial(i,j)+f_v(i,j)),
     &	density_liquid(i,j),v_l_r(kk),particle_diameter(i,j),
     &	dynamic_viscosity(i,j),velocity_circumferential(i,j),
     &	cylindrical_r(i,j),omiga(i,j))*dvdtr(density_vapor(i,j),
     &	(velocity_radial(i,j)+f_v(i,j)),density_liquid(i,j),v_l_r(kk),
     &	particle_diameter(i,j),dynamic_viscosity(i,j),
     &	velocity_circumferential(i,j),cylindrical_r(i,j),omiga(i,j))
	k2=delta_r(density_vapor(i,j),(velocity_radial(i,j)+f_v(i,j)),
     &	density_liquid(i,j),(v_l_r(kk)+0.5*k1),particle_diameter(i,j),
     &	dynamic_viscosity(i,j),velocity_circumferential(i,j),
     &	cylindrical_r(i,j),omiga(i,j))*dvdtr(density_vapor(i,j),
     &	(velocity_radial(i,j)+f_v(i,j)),density_liquid(i,j),
     &	(v_l_r(kk)+0.5*k1),particle_diameter(i,j),
     &	dynamic_viscosity(i,j),velocity_circumferential(i,j),
     &	cylindrical_r(i,j),omiga(i,j))
	k3=delta_r(density_vapor(i,j),(velocity_radial(i,j)+f_v(i,j)),
     &	density_liquid(i,j),(v_l_r(kk)+0.5*k2),particle_diameter(i,j),
     &	dynamic_viscosity(i,j),velocity_circumferential(i,j),
     &	cylindrical_r(i,j),omiga(i,j))*dvdtr(density_vapor(i,j),
     &	(velocity_radial(i,j)+f_v(i,j)),density_liquid(i,j),
     &	(v_l_r(kk)+0.5*k2),particle_diameter(i,j),
     &	dynamic_viscosity(i,j),velocity_circumferential(i,j),
     &	cylindrical_r(i,j),omiga(i,j))
	k4=delta_r(density_vapor(i,j),(velocity_radial(i,j)+f_v(i,j)),
     &	density_liquid(i,j),(v_l_r(kk)+k3),particle_diameter(i,j),
     &	dynamic_viscosity(i,j),velocity_circumferential(i,j),
     &	cylindrical_r(i,j),omiga(i,j))*dvdtr(density_vapor(i,j),
     &	(velocity_radial(i,j)+f_v(i,j)),density_liquid(i,j),
     &	(v_l_r(kk)+k3),particle_diameter(i,j),dynamic_viscosity(i,j),
     &	velocity_circumferential(i,j),cylindrical_r(i,j),omiga(i,j))
	v_l_r(kk+1)=v_l_r(kk)+(k1+2.0*k2+2.0*k3+k4)/6.0
	kk=kk+1
270	continue
	v_d_r(i,j)=v_l_r(kk)
250	continue
300	continue
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	��������һ��Һ���ٶȳ�
!	��һ��Һ�εľ����ٶ���Ϊ�����ٶȵ�ԭʼ����
	do 200 j=0,NSpan-1
	do 180 i=bd(j)+1,2*NStage*NStream-1
	v_l_c(0)=velocity_circumferential(i,j)
	v_l_c(1)=velocity_circumferential(i,j)+0.1
!	�����ʼ�ٶ�
	kk=1
	do 150 while (abs(v_l_c(kk)-v_l_c(kk-1))>1D-6)
	k1=delta_c(density_vapor(i,j),
     &	(velocity_circumferential(i,j)+f_v(i,j)),density_liquid(i,j),
     &	v_l_c(kk),particle_diameter(i,j),dynamic_viscosity(i,j),
     &	v_d_r(i,j),cylindrical_r(i,j),omiga(i,j))*
     &	dvdtc(density_vapor(i,j),
     &	(velocity_circumferential(i,j)+f_v(i,j)),density_liquid(i,j),
     &	v_l_c(kk),particle_diameter(i,j),dynamic_viscosity(i,j),
     &	v_d_r(i,j),cylindrical_r(i,j),omiga(i,j))
	k2=delta_c(density_vapor(i,j),
     &	(velocity_circumferential(i,j)+f_v(i,j)),density_liquid(i,j),
     &	(v_l_c(kk)+0.5*k1),particle_diameter(i,j),
     &	dynamic_viscosity(i,j),v_d_r(i,j),cylindrical_r(i,j),
     &	omiga(i,j))*dvdtc(density_vapor(i,j),
     &	(velocity_circumferential(i,j)+f_v(i,j)),density_liquid(i,j),
     &	v_l_c(kk),particle_diameter(i,j),dynamic_viscosity(i,j),
     &	v_d_r(i,j),cylindrical_r(i,j),omiga(i,j))
	k3=delta_c(density_vapor(i,j),
     &	(velocity_circumferential(i,j)+f_v(i,j)),density_liquid(i,j),
     &	(v_l_c(kk)+0.5*k2),particle_diameter(i,j),
     &	dynamic_viscosity(i,j),v_d_r(i,j),cylindrical_r(i,j),
     &	omiga(i,j))*dvdtc(density_vapor(i,j),
     &	(velocity_circumferential(i,j)+f_v(i,j)),density_liquid(i,j),
     &	v_l_c(kk),particle_diameter(i,j),dynamic_viscosity(i,j),
     &	v_d_r(i,j),cylindrical_r(i,j),omiga(i,j))
	k4=delta_c(density_vapor(i,j),
     &	(velocity_circumferential(i,j)+f_v(i,j)),density_liquid(i,j),
     &	(v_l_c(kk)+k3),particle_diameter(i,j),dynamic_viscosity(i,j),
     &	v_d_r(i,j),cylindrical_r(i,j),omiga(i,j))*
     &	dvdtc(density_vapor(i,j),
     &	(velocity_circumferential(i,j)+f_v(i,j)),density_liquid(i,j),
     &	v_l_c(kk),particle_diameter(i,j),dynamic_viscosity(i,j),
     &	v_d_r(i,j),cylindrical_r(i,j),omiga(i,j))
	v_l_c(kk+1)=v_l_c(kk)+(k1+2.0*k2+2.0*k3+k4)/6.0
	kk=kk+1
150	continue
	v_d_c(i,j)=v_l_c(kk)
180	continue
200	continue
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	����һ��Һ�������ٶȳ�
	do 100, k=0,NSpan-1
!	����NSpan��Kʱ��������ٶȳ�
	istart=bd(k)
!	���ڿ�ʼ����
	v_l_a(0)=velocity_axial(istart,k)
	p_l_a(0)=cylindrical_z(istart,k)
!	�����ʼ�ٶȺ�λ��
	i=1
	do 40 while (p_l_a(i-1)<=
     &	cylindrical_z(NStream*(istart/NStream+1)-1,k))
	do j=0,NStream-1
		qx_weizhi(j)=cylindrical_z(NStream*(istart/NStream)+j,k)
		qx_sudu(j)=velocity_axial(NStream*(istart/NStream)+j,k)+
     &		f_v(NStream*(istart/NStream)+j,k)
		qx_midu(j)=density_vapor(NStream*(istart/NStream)+j,k)
		qx_nianxing(j)=dynamic_viscosity(NStream*(istart/NStream)+j,k)
		kl_midu(j)=density_liquid(NStream*(istart/NStream)+j,k)
		kl_zhijing(j)=particle_diameter(NStream*(istart/NStream)+j,k)
	enddo
	call SPLINE(qx_weizhi,qx_sudu,NStream-1,p_l_a(i-1),
     &	qxsudu,temp)
	call SPLINE(qx_weizhi,qx_midu,NStream-1,p_l_a(i-1),
     &	qxmidu,temp)
	call SPLINE(qx_weizhi,qx_nianxing,NStream-1,p_l_a(i-1),
     &	qxnianxing,temp)
	call SPLINE(qx_weizhi,kl_midu,NStream-1,p_l_a(i-1),
     &	klmidu,temp)
	call SPLINE(qx_weizhi,kl_zhijing,NStream-1,p_l_a(i-1),
     &	klzhijing,temp)
	k1=delta_z(qxmidu,qxsudu,klmidu,v_l_a(i-1),klzhijing,
     &	qxnianxing)*dvdtz(qxmidu,qxsudu,klmidu,
     &	v_l_a(i-1),klzhijing,qxnianxing)
	k2=delta_z(qxmidu,qxsudu,klmidu,v_l_a(i-1),klzhijing,
     &	qxnianxing)*dvdtz(qxmidu,qxsudu,klmidu,
     &	(v_l_a(i-1)+0.5*k1),klzhijing,qxnianxing)
	k3=delta_z(qxmidu,qxsudu,klmidu,v_l_a(i-1),klzhijing,
     &	qxnianxing)*dvdtz(qxmidu,qxsudu,klmidu,
     &	(v_l_a(i-1)+0.5*k2),klzhijing,qxnianxing)
	k4=delta_z(qxmidu,qxsudu,klmidu,v_l_a(i-1),klzhijing,
     &	qxnianxing)*dvdtz(qxmidu,qxsudu,klmidu,
     &	(v_l_a(i-1)+k3),klzhijing,qxnianxing)
	v_l_a(i)=v_l_a(i-1)+(k1+2.0*k2+2.0*k3+k4)/6.0
	p_l_a(i)=p_l_a(i-1)+
     &	0.5*delta_z(qxmidu,qxsudu,klmidu,v_l_a(i-1),klzhijing,
     &	qxnianxing)*(v_l_a(i)+v_l_a(i-1))
	i=i+1
40	continue
	nend=i-1
	do i=istart-NStream*(istart/NStream),NStream-1
		do j=0,nend-1
	if (qx_weizhi(i)>=p_l_a(j).and.qx_weizhi(i)<=p_l_a(j+1)) then
		v_d_a(NStream*(istart/NStream)+i,k)=v_l_a(j)+
     &	(v_l_a(j+1)-v_l_a(j))*(qx_weizhi(i)-p_l_a(j))/
     &	(p_l_a(j+1)-p_l_a(j))
	endif
		enddo
	enddo
!	����NSpan��Kʱ���θ�������ٶȳ�
	if(istart/NStream+2<=2*NStage) then
	do 90, ii=istart/NStream+2,2*NStage
	v_l_a(0)=v_d_a((ii-1)*NStream-1,k)
	p_l_a(0)=cylindrical_z((ii-1)*NStream,k)
!	���θ�����ĳ���
	i=1
	do 80 while (p_l_a(i-1)<=
     &	cylindrical_z(ii*NStream-1,k))
	do j=0,NStream-1
		qx_weizhi(j)=cylindrical_z((ii-1)*NStream+j,k)
		qx_sudu(j)=velocity_axial((ii-1)*NStream+j,k)+
     &	f_v((ii-1)*NStream+j,k)
		qx_midu(j)=density_vapor((ii-1)*NStream+j,k)
		qx_nianxing(j)=dynamic_viscosity((ii-1)*NStream+j,k)
		kl_midu(j)=density_liquid((ii-1)*NStream+j,k)
		kl_zhijing(j)=particle_diameter((ii-1)*NStream+j,k)
	enddo
	call SPLINE(qx_weizhi,qx_sudu,NStream-1,p_l_a(i-1),
     &	qxsudu,temp)
	call SPLINE(qx_weizhi,qx_midu,NStream-1,p_l_a(i-1),
     &	qxmidu,temp)
	call SPLINE(qx_weizhi,qx_nianxing,NStream-1,p_l_a(i-1),
     &	qxnianxing,temp)
	call SPLINE(qx_weizhi,kl_midu,NStream-1,p_l_a(i-1),
     &	klmidu,temp)
	call SPLINE(qx_weizhi,kl_zhijing,NStream-1,p_l_a(i-1),
     &	klzhijing,temp)
	k1=delta_z(qxmidu,qxsudu,klmidu,v_l_a(i-1),klzhijing,
     &	qxnianxing)*dvdtz(qxmidu,qxsudu,klmidu,
     &	v_l_a(i-1),klzhijing,qxnianxing)
	k2=delta_z(qxmidu,qxsudu,klmidu,v_l_a(i-1),klzhijing,
     &	qxnianxing)*dvdtz(qxmidu,qxsudu,klmidu,
     &	(v_l_a(i-1)+0.5*k1),klzhijing,qxnianxing)
	k3=delta_z(qxmidu,qxsudu,klmidu,v_l_a(i-1),klzhijing,
     &	qxnianxing)*dvdtz(qxmidu,qxsudu,klmidu,
     &	(v_l_a(i-1)+0.5*k2),klzhijing,qxnianxing)
	k4=delta_z(qxmidu,qxsudu,klmidu,v_l_a(i-1),klzhijing,
     &	qxnianxing)*dvdtz(qxmidu,qxsudu,klmidu,
     &	(v_l_a(i-1)+k3),klzhijing,qxnianxing)
	v_l_a(i)=v_l_a(i-1)+(k1+2.0*k2+2.0*k3+k4)/6.0
	p_l_a(i)=p_l_a(i-1)+
     &	0.5*delta_z(qxmidu,qxsudu,klmidu,v_l_a(i-1),klzhijing,
     &	qxnianxing)*(v_l_a(i)+v_l_a(i-1))
	i=i+1
80	continue
	nend=i-1
	do i=0,NStream-1
		do j=0,nend
		if(qx_weizhi(i)>=p_l_a(j).and.qx_weizhi(i)<p_l_a(j+1))then
	v_d_a((ii-1)*NStream+i,k)=v_l_a(j)+(v_l_a(j+1)-v_l_a(j))*
     &	(qx_weizhi(i)-p_l_a(j))/(p_l_a(j+1)-p_l_a(j))
		endif
		enddo
	enddo
90	continue
	else
	endif
100	continue
!	һ��Һ�������ٶȳ�������ϣ�������v_d_a(i,k)��
	print *,'һ��Һ���ٶȼ������!'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	�������Һ�������ٶȳ�
	print *,'��ʼ�������Һ���ٶ�......'
	call stator_rotor_a_r(slr,slz,str,stz,rlr,rlz,rtr,rtz)
!	�����Ӻ������������Ҷ�Ͷ�Ҷǰ��Ե����
	do i=1,NStage
		do j=0,NSpan-1
			chushisudu(i,j)=0.3D0
			chushiweizhi(i,j)=stz(i,j)
		enddo
	enddo
	do 101, k=0,NSpan-1
!	����NSpan��Kʱ��������ٶȳ�
	do 50, ii=1,NStage
	v_l_a(0)=chushisudu(ii,k)
	p_l_a(0)=chushiweizhi(ii,k)
!	�������ڳ�ʼ�ٶȺͳ�ʼλ��
	ibjcssy=2*(ii-1)*NStream
	bjzlwz=cylindrical_z((2*ii-1)*NStream-1,k)
!	ibjcssy:���ڳ�ʼ����,bjzlwz:��������λ��
	i=1
	do 41 while (p_l_a(i-1)<=bjzlwz)
	do j=0,NStream-1
		qx_weizhi(j)=cylindrical_z(ibjcssy+j,k)
		qx_sudu(j)=velocity_axial(ibjcssy+j,k)
		qx_midu(j)=density_vapor(ibjcssy+j,k)
		qx_nianxing(j)=dynamic_viscosity(ibjcssy+j,k)
		kl_midu(j)=density_liquid(ibjcssy+j,k)
		kl_surfacetension(j)=surface_tension(ibjcssy+j,k)
	enddo
	call SPLINE(qx_weizhi,qx_midu,NStream-1,p_l_a(0),
     &	chushiqxmidu,temp)
	call SPLINE(qx_weizhi,qx_sudu,NStream-1,p_l_a(0),
     &	chushiqxsudu,temp)
	call SPLINE(qx_weizhi,kl_surfacetension,NStream-1,p_l_a(0),
     &	klsurfacetension,temp)
	call SPLINE(qx_weizhi,kl_midu,NStream-1,p_l_a(0),
     &	kldmidu(ii,k),temp)
!	�����ܶȺͿ������������ɳ�ʼ����ȷ�����ڱ����ļ����б��ֲ���
    	kldzhijing(ii,k)=20.0*klsurfacetension/chushiqxmidu/
     &	((chushiqxsudu-chushisudu(ii,k))**2)/2.0
!	�����ٽ�Τ����ȷ�����ο���ֱ������ȡ�ٽ�ֱ����һ����Ϊƽ������ֱ��
	call SPLINE(qx_weizhi,qx_sudu,NStream-1,p_l_a(i-1),
     &	qxsudu,temp)
	call SPLINE(qx_weizhi,qx_midu,NStream-1,p_l_a(i-1),
     &	qxmidu,temp)
	call SPLINE(qx_weizhi,qx_nianxing,NStream-1,p_l_a(i-1),
     &	qxnianxing,temp)
	k1=delta_z(qxmidu,qxsudu,kldmidu(ii,k),v_l_a(i-1),
     &	kldzhijing(ii,k),qxnianxing)*dvdtz(qxmidu,qxsudu,
     &	kldmidu(ii,k),v_l_a(i-1),kldzhijing(ii,k),
     &	qxnianxing)/1000.0
	k2=delta_z(qxmidu,qxsudu,kldmidu(ii,k),v_l_a(i-1),
     &	kldzhijing(ii,k),qxnianxing)*dvdtz(qxmidu,qxsudu,
     &	kldmidu(ii,k),(v_l_a(i-1)+0.5*k1),kldzhijing(ii,k),
     &	qxnianxing)/1000.0
	k3=delta_z(qxmidu,qxsudu,kldmidu(ii,k),v_l_a(i-1),
     &	kldzhijing(ii,k),qxnianxing)*dvdtz(qxmidu,qxsudu,
     &	kldmidu(ii,k),(v_l_a(i-1)+0.5*k2),kldzhijing(ii,k),
     &	qxnianxing)/1000.0
	k4=delta_z(qxmidu,qxsudu,kldmidu(ii,k),v_l_a(i-1),
     &	kldzhijing(ii,k),qxnianxing)*dvdtz(qxmidu,qxsudu,
     &	kldmidu(ii,k),(v_l_a(i-1)+k3),kldzhijing(ii,k),
     &	qxnianxing)/1000.0
	v_l_a(i)=v_l_a(i-1)+(k1+2.0*k2+2.0*k3+k4)/6.0
	p_l_a(i)=p_l_a(i-1)+0.5*delta_z(qxmidu,qxsudu,kldmidu(ii,k),
     &	v_l_a(i-1),kldzhijing(ii,k),qxnianxing)*(v_l_a(i)+
     &	v_l_a(i-1))/1000.0
	i=i+1
41	continue
	nend=i-1
	do i=0,NStream-1
		do j=0,nend
		if(qx_weizhi(i)>=p_l_a(j).and.qx_weizhi(i)<p_l_a(j+1))then
		v_ds_a(ibjcssy+i,k)=v_l_a(j)+(v_l_a(j+1)-v_l_a(j))*
     &	(qx_weizhi(i)-p_l_a(j))/(p_l_a(j+1)-p_l_a(j))
		endif
		enddo
	enddo
!	��Ҷ������������ٶȼ������
!	���㶯Ҷ������������ٶ�
!		ibjcssy=���ڳ�ʼ����
	v_l_a(0)=v_ds_a(ibjcssy+NStream-1,k)
	p_l_a(0)=cylindrical_z(ibjcssy+NStream,k)
!	��Ҷ����������ٶȺ;�Ҷ������ڵ�һ��
	i=1
	do 81 while (p_l_a(i-1)<=
     &	cylindrical_z(ibjcssy+2*NStream-1,k))
	do j=0,NStream-1
		qx_weizhi(j)=cylindrical_z(ibjcssy+NStream+j,k)
		qx_sudu(j)=velocity_axial(ibjcssy+NStream+j,k)
		qx_midu(j)=density_vapor(ibjcssy+NStream+j,k)
		qx_nianxing(j)=dynamic_viscosity(ibjcssy+NStream+j,k)
		kl_midu(j)=density_liquid(ibjcssy+NStream,k)
		kl_surfacetension(j)=surface_tension(ibjcssy+NStream+j,k)
	enddo
	call SPLINE(qx_weizhi,qx_sudu,NStream-1,p_l_a(i-1),
     &	qxsudu,temp)
	call SPLINE(qx_weizhi,qx_midu,NStream-1,p_l_a(i-1),
     &	qxmidu,temp)
	call SPLINE(qx_weizhi,qx_nianxing,NStream-1,p_l_a(i-1),
     &	qxnianxing,temp)
	k1=delta_z(qxmidu,qxsudu,kldmidu(ii,k),v_l_a(i-1),
     &	kldzhijing(ii,k),qxnianxing)*dvdtz(qxmidu,qxsudu,
     &	kldmidu(ii,k),v_l_a(i-1),kldzhijing(ii,k),
     &	qxnianxing)/1000.0
	k2=delta_z(qxmidu,qxsudu,kldmidu(ii,k),v_l_a(i-1),
     &	kldzhijing(ii,k),qxnianxing)*dvdtz(qxmidu,qxsudu,
     &	kldmidu(ii,k),(v_l_a(i-1)+0.5*k1),kldzhijing(ii,k),
     &	qxnianxing)/1000.0
	k3=delta_z(qxmidu,qxsudu,kldmidu(ii,k),v_l_a(i-1),
     &	kldzhijing(ii,k),qxnianxing)*dvdtz(qxmidu,qxsudu,
     &	kldmidu(ii,k),(v_l_a(i-1)+0.5*k2),kldzhijing(ii,k),
     &	qxnianxing)/1000.0
	k4=delta_z(qxmidu,qxsudu,kldmidu(ii,k),v_l_a(i-1),
     &	kldzhijing(ii,k),qxnianxing)*dvdtz(qxmidu,qxsudu,
     &	kldmidu(ii,k),(v_l_a(i-1)+k3),kldzhijing(ii,k),
     &	qxnianxing)/1000.0
	v_l_a(i)=v_l_a(i-1)+(k1+2.0*k2+2.0*k3+k4)/6.0
	p_l_a(i)=p_l_a(i-1)+0.5*delta_z(qxmidu,qxsudu,kldmidu(ii,k),
     &	v_l_a(i-1),kldzhijing(ii,k),qxnianxing)*(v_l_a(i)+
     &	v_l_a(i-1))/1000.0
	i=i+1
81	continue
	nend=i-1
	do i=0,NStream-1
		do j=0,nend
		if(qx_weizhi(i)>=p_l_a(j).and.qx_weizhi(i)<p_l_a(j+1))then
	v_ds_a(ibjcssy+NStream+i,k)=v_l_a(j)+(v_l_a(j+1)-v_l_a(j))*
     &	(qx_weizhi(i)-p_l_a(j))/(p_l_a(j+1)-p_l_a(j))
		endif
		enddo
	enddo
50	continue
101	continue
!	����Һ�������ٶȳ�������ϣ�������v_ds_a(i,k)��
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	�������Һ�������ٶȳ�
!	�ھ�Ҷ�Ͷ�Ҷ֮���Blade-to-Blade��������,
!	ֻ�й�������������������,���������ƶ������Ӱ����Ժ���,
!	��˼�������������ڶ���Һ���ٶȵķ���������ٶȵķ���һ��,
!	�����Ա���,���Լ򻯶���Һ�������ٶȵļ���.
	do i=0,2*NStage*NStream-1
	do j=0,NSpan-1
		v_ds_c(i,j)=v_ds_a(i,j)*
     &	velocity_circumferential(i,j)/velocity_axial(i,j)
	enddo
	enddo
!	����Һ�������ٶȳ�������ϣ�������v_ds_c(i,k)��
	print *,'�������Һ���ٶȼ������!'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	Ϊ�˼��㷽�㣬�ȼ���һ��Һ���ٶȳ����������֮�����滻���ڼ���
c	do i=0,2*NStage*NStream-1
c	do j=0,NSpan-1
c		v_d_r(i,j)=velocity_radial(i,j)+0.01
c		v_d_c(i,j)=velocity_circumferential(i,j)-0.01
c		v_d_a(i,j)=velocity_axial(i,j)-0.2
c	enddo
c	enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	��������ѧ��ʧ
	call thermodynamicloss()
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	������μд���ʧ
	call reentainmentloss()
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	�����ƶ���ʧ
	call impactloss()
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	������ˮ��ʧ
	call collectedwaterloss()
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	����������ʧ
	call centrifugingloss()
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	����Ҷդ�����ʺ���
!	��ȷ
!	��֤����Ҷդ��Ҷ��Һ������
!	��ȷ
c	call DEPOSITION_RATE(dp_rt_s,dp_rt_r)
c	call FLOWRATE(liquid_s_massflow,liquid_r_massflow,
c     &	vapor_s_massflow,vapor_r_massflow)
c	do k=1,NStage
c		do j=0,NSpan-1
c			print *,j,dp_rt_s(k,j),dp_rt_r(k,j)
c	print *,k,j,vapor_s_massflow(k,j),vapor_r_massflow(k,j)
c		enddo
c	enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	end
