!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	subroutine data_order(a,b)
	parameter (NStage=10,NSpan=11,NStream=11)
	double precision a(2*NStage),tmp(2*NStage)
	integer ip,b(2*NStage),bignumber
	bignumber=2*NStage*NSpan*NStream
	do i=1,2*NStage
		tmp(i)=a(i)
	enddo
	do j=1,2*NStage
		ip=1
		do i=1,2*NStage
			if (a(i) < a(ip)) then
				ip=i
			endif
		enddo
		a(ip)=bignumber
		b(j)=ip
	enddo
	do i=1,2*NStage
		a(i)=tmp(i)
	enddo
	end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	��������������С�����˳�������������
	subroutine adjust_data_order(NStage,NSpan,NStream,NNode,z,zb,za)
	integer NStage,NSpan,NStream,NNode
	integer NI(2*NStage),bbb(2*NStage)
	double precision tmp(2*NStage)
	double precision z(0:NNode-1),zb(0:NNode-1),za(0:NNode-1)
	do k=1,2*NStage
		NI(k)=NSpan*NStream*(k-1)
		tmp(k)=z(NI(k))
	enddo
	call data_order(tmp,bbb)
	do k=1,2*NStage
		do i=0,NSpan*NStream-1
			j=(k-1)*NSpan*NStream+i
			za(j)=zb((bbb(k)-1)*NSpan*NStream+i)
		enddo
	enddo
	end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	��ֱ������ϵת��ΪԲ������ϵ
	subroutine adjust_coordinate(NStage,NSpan,NStream,NNode,z,zz)
	integer NStage,NSpan,NStream,NNode
	double precision z(0:NNode-1),zz(0:2*NStage*NStream-1,0:NSpan-1)
	do i=0,2*NStage*NStream-1
		do j=0,NSpan-1
			zz(i,j)=z(NSpan*i+j)
		enddo
	enddo
	end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	�������������������ֵ��һ�ε���
	SUBROUTINE SPLINE(X,Y,N,SX,SY,DY)
      double precision X(0:N),Y(0:N),H(N),SX,SY,DY
	double precision LAMDA(N-1),MU(N-1),D(N-1),M(0:N)
 		IF (N<1) THEN
	STOP 'ERROR! THE NUMBER OF ARRAY IS WRONG!'
		ELSEIF (N==1) THEN
	SY=(SX-X(1))/(X(0)-X(1))*Y(0)+(SX-X(0))/(X(1)-X(0))*Y(1)
	DY=(Y(1)-Y(0))/(X(1)-X(0))
		ELSEIF (N==2) THEN
	SY= (SX-X(1))*(SX-X(2))/(X(0)-X(1))/(X(0)-X(2))*Y(0)+
     &	(SX-X(0))*(SX-X(2))/(X(1)-X(0))/(X(1)-X(2))*Y(1)+
     &	(SX-X(0))*(SX-X(1))/(X(2)-X(0))/(X(2)-X(1))*Y(2)
	DY= (SX-X(1)+SX-X(2))/(X(0)-X(1))/(X(0)-X(2))*Y(0)+
     &	(SX-X(0)+SX-X(2))/(X(1)-X(0))/(X(1)-X(2))*Y(1)+
     &	(SX-X(0)+SX-X(1))/(X(2)-X(0))/(X(2)-X(1))*Y(2)
		ELSE
      DO J=1,N
		H(J)=X(J)-X(J-1)
c		IF (H(J)<=1D-8)	H(J)=1D-8
!	����������ĳ�������ݼ����С���غϣ�����ำ��Сֵ����֤�����в��������
	ENDDO
      DO J=1,N-1
		LAMDA(J)=H(J+1)/(H(J)+H(J+1))
		MU(J)=1.0-LAMDA(J)
		D(J)=(((Y(J+1)-Y(J))/H(J+1))-((Y(J)-Y(J-1))/H(J)))
     &	*6/(H(J+1)+H(J))
      ENDDO
      CALL TDMA(MU,LAMDA,D,N,M)
	IF(SX.LT.X(0)) THEN
		SY=(M(0)*(X(1)-SX)**3+M(1)*(SX-X(0))**3)/6/H(1)
		SY=SY+(Y(0)-M(0)*H(1)**2/6.0)*(X(1)-SX)/H(1)
		SY=SY+(Y(1)-M(1)*H(1)**2/6.0)*(SX-X(0))/H(1)
		DY=(-((X(1)-SX)**2)*M(0)+((SX-X(0))**2)*M(1))/H(1)/2.0
		DY=DY+(Y(1)-Y(0))/H(1)+(M(0)-M(1))*H(1)/6.0
	ELSE IF(SX.GE.X(N)) THEN
		SY=(M(N-1)*(X(N)-SX)**3+M(N)*(SX-X(N-1))**3)/6/H(N)
		SY=SY+(Y(N-1)-M(N-1)*H(N)**2/6.0)*(X(N)-SX)/H(N)
		SY=SY+(Y(N)-M(N)*H(N)**2/6.0)*(SX-X(N-1))/H(N)
		DY=(-((X(N)-SX)**2)*M(N-1)+((SX-X(N-1))**2)*M(N))/H(N)/2.0
		DY=DY+(Y(N)-Y(N-1))/H(N)+(M(N-1)-M(N))*H(N)/6.0	
	ELSE 
		DO J=1,N
			IF ((SX.GE.X(J-1)).AND.(SX.LT.X(J))) THEN
			SY=(M(J-1)*(X(J)-SX)**3+M(J)*(SX-X(J-1))**3)/6/H(J)
			SY=SY+(Y(J-1)-M(J-1)*H(J)**2/6.0)*(X(J)-SX)/H(J)
			SY=SY+(Y(J)-M(J)*H(J)**2/6.0)*(SX-X(J-1))/H(J)
			DY=(-((X(J)-SX)**2)*M(J-1)+((SX-X(J-1))**2)*M(J))/H(J)/2.0
			DY=DY+(Y(J)-Y(J-1))/H(J)+(M(J-1)-M(J))*H(J)/6.0
			ENDIF
		ENDDO
      END IF
		ENDIF
      END
      SUBROUTINE TDMA(MU,LAMDA,D,N,M)
      double precision MU(N-1),LAMDA(N-1),D(N-1),M(0:N)
	double precision A(2:N-1),B(N-1),C(N-2)
	double precision L(N-1),R(N-2),X(N-1),Y(N-1)
	M(0)=0.0
	M(N)=0.0
!	������Ȼ���������������壬�˵㴦�Ķ��ε���Ϊ0��
	DO I=1,N-2
		A(I+1)=MU(I+1)
		C(I)=LAMDA(I)
	ENDDO
	DO I=1,N-1
		B(I)=2.0
	ENDDO
	L(1)=B(1)
	Y(1)=D(1)/L(1)
	DO I=2,N-1
		R(I-1)=C(I-1)/L(I-1)
		L(I)=B(I)-A(I)*R(I-1)
		Y(I)=(D(I)-A(I)*Y(I-1))/L(I)
	ENDDO
	X(N-1)=Y(N-1)
	DO I=N-2,1,-1
		X(I)=Y(I)-R(I)*X(I+1)
	ENDDO
	DO I=1,N-1
		M(I)=X(I)
	ENDDO
	END
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	�����������������ֱ�߽�������
!	������y1>=y2�������½��������ȷ��
	SUBROUTINE SPLINE_INTERSECTION(x,y,n,x1,y1,x2,y2,x3,y3)
	integer n
	double precision x(0:n),y(0:n)
	double precision x1,y1,x2,y2,x3,y3
	double precision beta
	double precision xx(0:n),yy(0:n)
	double precision x1n,y1n,x2n,y2n
	double precision x3temp,y3temp,temp
	beta=asin((x1-x2)/sqrt((x1-x2)**2+(y1-y2)**2))
	do i=0,n
		xx(i)=(x(i)-x(0))*cos(beta)-(y(i)-y(0))*sin(beta)+x(0)
		yy(i)=(x(i)-x(0))*sin(beta)+(y(i)-y(0))*cos(beta)+y(0)
	enddo
	x1n=(x1-x(0))*cos(beta)-(y1-y(0))*sin(beta)+x(0)
	y1n=(x1-x(0))*sin(beta)+(y1-y(0))*cos(beta)+y(0)
	x2n=(x2-x(0))*cos(beta)-(y2-y(0))*sin(beta)+x(0)
	y2n=(x2-x(0))*sin(beta)+(y2-y(0))*cos(beta)+y(0)
	x3temp=x1n
	x3temp=x2n
	call SPLINE(xx,yy,n,x3temp,y3temp,temp)
	x3=(x3temp-x(0))*cos(-beta)-(y3temp-y(0))*sin(-beta)+x(0)
	y3=(x3temp-x(0))*sin(-beta)+(y3temp-y(0))*cos(-beta)+y(0)
	END
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	�Լ�����ı߽����ݽ��е���,�߽紦�����������
	SUBROUTINE DATA_ADJUST_COORDINATE(NStage,NSpan,NStream,xyz)
	double precision xyz(0:2*NStage*NStream-1,0:NSpan-1)
	double precision temp,eps
	eps=1D-6
	do i=1,2*NStage-1
		do j=0,NSpan-1
			temp=(xyz(i*NStream-1,j)+xyz(i*NStream,j))*0.5
			xyz(i*NStream-1,j)=temp-eps
			xyz(i*NStream,j)=temp+eps
		enddo
	enddo
	END
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	�Լ�����ı߽����ݽ��е���,�߽紦��ֵ��ͬ
	SUBROUTINE DATA_ADJUST_VALUE(NStage,NSpan,NStream,xyz)
	double precision xyz(0:2*NStage*NStream-1,0:NSpan-1)
	double precision temp
	do i=1,2*NStage-1
		do j=0,NSpan-1
			temp=(xyz(i*NStream-1,j)+xyz(i*NStream,j))*0.5
			xyz(i*NStream-1,j)=temp
			xyz(i*NStream,j)=temp
		enddo
	enddo
	END
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	���ж�Ԫ���ݲ�ֵ
	SUBROUTINE SPLINE_AREA(zfunction,x0,y0,z0)
	implicit double precision (a-h,o-z)
	parameter (NStage=10,NSpan=11,NStream=11)
	common/com01/cylindrical_r(0:2*NStage*NStream-1,0:NSpan-1)
	common/com02/cylindrical_z(0:2*NStage*NStream-1,0:NSpan-1)
	double precision zfunction(0:2*NStage*NStream-1,0:NSpan-1)
	double precision x0,y0,z0
	double precision xx(0:2*NStage*NStream-1)
	double precision yy(0:2*NStage*NStream-1)
	double precision zz(0:2*NStage*NStream-1)
	double precision yyy(0:NSpan-1)
	double precision zzz(0:NSpan-1)
	double precision temp
	do j=0,NSpan-1
		do i=0,2*NStage*NStream-1
			xx(i)=cylindrical_z(i,j)
			yy(i)=cylindrical_r(i,j)
			zz(i)=zfunction(i,j)
		enddo
		call SPLINE(xx,yy,2*NStage*NStream-1,x0,yyy(j),temp)
		call SPLINE(xx,zz,2*NStage*NStream-1,x0,zzz(j),temp)
	enddo
	call SPLINE(yyy,zzz,NSpan-1,y0,z0,temp)
	END
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	����Һ�������ٶ�ƫ����
	function dvdtz(a,b,c,d,e,f)
	double precision a,b,c,d,e,f
	double precision re,cd,qx_midu,qx_sudu,kl_midu,kl_sudu,kl_zhijing,
     &	qx_nianxing,dvdtz
	qx_midu=a
	qx_sudu=b
	kl_midu=c
	kl_sudu=d
	kl_zhijing=e
	qx_nianxing=f
	re=qx_midu*abs(qx_sudu-kl_sudu)*kl_zhijing/qx_nianxing
	if (re<1D-8) then
		dvdtz=18*qx_nianxing*(qx_sudu-kl_sudu)/
     &		kl_zhijing/kl_zhijing/kl_midu
	elseif (re<0.1) then
		cd=24.0/re
		dvdtz=3*cd*qx_midu*(qx_sudu-kl_sudu)*abs(qx_sudu-kl_sudu)
     &	/4/kl_zhijing/kl_midu
	elseif(re<=1.0) then
		cd=22.73/re+0.0903/re/re+3.69
		dvdtz=3*cd*qx_midu*(qx_sudu-kl_sudu)*abs(qx_sudu-kl_sudu)
     &	/4/kl_zhijing/kl_midu
	elseif(re<=10.0) then
		cd=29.1667/re-3.8889/re/re+1.222
		dvdtz=3*cd*qx_midu*(qx_sudu-kl_sudu)*abs(qx_sudu-kl_sudu)
     &	/4/kl_zhijing/kl_midu
	elseif(re<=100.0) then
		cd=46.5/re-116.67/re/re+0.6167
		dvdtz=3*cd*qx_midu*(qx_sudu-kl_sudu)*abs(qx_sudu-kl_sudu)
     &	/4/kl_zhijing/kl_midu
	elseif(re<=1000.0) then
		cd=98.33/re-2778/re/re+0.3644
		dvdtz=3*cd*qx_midu*(qx_sudu-kl_sudu)*abs(qx_sudu-kl_sudu)
     &	/4/kl_zhijing/kl_midu
	elseif(re<=5000.0) then
		cd=148.62/re-47500/re/re+0.357
		dvdtz=3*cd*qx_midu*(qx_sudu-kl_sudu)*abs(qx_sudu-kl_sudu)
     &	/4/kl_zhijing/kl_midu
	elseif(re<=10000.0) then
		cd=-490.546/re+578700/re/re+0.46
		dvdtz=3*cd*qx_midu*(qx_sudu-kl_sudu)*abs(qx_sudu-kl_sudu)
     &	/4/kl_zhijing/kl_midu	
	elseif(re<=50000.0) then
		cd=-1662.5/re+5416700/re/re+0.5191
		dvdtz=3*cd*qx_midu*(qx_sudu-kl_sudu)*abs(qx_sudu-kl_sudu)
     &	/4/kl_zhijing/kl_midu
	endif
	end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	����Һ�����򵱵�ʱ�䲽��
	function delta_z(a,b,c,d,e,f)
	double precision a,b,c,d,e,f
	double precision qx_midu,qx_sudu,kl_midu,kl_sudu,kl_zhijing,qx_nianxing
	double precision delta_z
	qx_midu=a
	qx_sudu=b
	kl_midu=c
	kl_sudu=d
	kl_zhijing=e
	qx_nianxing=f
	if (abs(qx_sudu-kl_sudu)<1D-8) then
		delta_z=kl_zhijing*kl_zhijing*kl_midu/18.0/qx_nianxing
	else
		delta_z=1.0/dvdtz(qx_midu,qx_sudu,kl_midu,kl_sudu,
     &	kl_zhijing,qx_nianxing)*(qx_sudu-kl_sudu)
	endif
!	��600MW���ֻ����������ϵ��Ϊ2.7
	delta_z=delta_z*2.0
c	delta_z=delta_z*2.5
!	�������ٶȺ�λ�ñ�������ʱ������Ӷ��ı䣬�������ʱ�䲽��̫С��С��ʵ�б����ľ���
!	������ϵͳֱ�Ӻ��ԣ���˿��Խ�ʱ�䲽���ʵ�����!
	end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	����Һ�������ٶ�ƫ����
	function dvdtc(a,b,c,d,e,f,cvr,cr,co)
	double precision a,b,c,d,e,f,cvr,cr,co
	double precision re,cd,qx_midu,qx_sudu,kl_midu,kl_sudu,kl_zhijing,
     &	qx_nianxing,kl_r_sudu,ban_jing,jiao_sudu,dvdtc
	qx_midu=a
	qx_sudu=b
	kl_midu=c
	kl_sudu=d
	kl_zhijing=e
	qx_nianxing=f
	kl_r_sudu=cvr
	ban_jing=cr
	jiao_sudu=co
	re=qx_midu*abs(qx_sudu-kl_sudu)*kl_zhijing/qx_nianxing
	if (re<1D-8) then
		dvdtc=18*qx_nianxing*(qx_sudu-kl_sudu)/kl_zhijing/
     &	kl_zhijing/kl_midu-2*kl_r_sudu*(kl_sudu/ban_jing+jiao_sudu)
	elseif (re<0.1) then
		cd=24.0/re
		dvdtc=3*cd*qx_midu*(qx_sudu-kl_sudu)*abs(qx_sudu-kl_sudu)
     &	/4/kl_zhijing/kl_midu-2*kl_r_sudu*(kl_sudu/ban_jing+jiao_sudu)
	elseif(re<=1.0) then
		cd=22.73/re+0.0903/re/re+3.69
		dvdtc=3*cd*qx_midu*(qx_sudu-kl_sudu)*abs(qx_sudu-kl_sudu)
     &	/4/kl_zhijing/kl_midu-2*kl_r_sudu*(kl_sudu/ban_jing+jiao_sudu)
	elseif(re<=10.0) then
		cd=29.1667/re-3.8889/re/re+1.222
		dvdtc=3*cd*qx_midu*(qx_sudu-kl_sudu)*abs(qx_sudu-kl_sudu)
     &	/4/kl_zhijing/kl_midu-2*kl_r_sudu*(kl_sudu/ban_jing+jiao_sudu)
	elseif(re<=100.0) then
		cd=46.5/re-116.67/re/re+0.6167
		dvdtc=3*cd*qx_midu*(qx_sudu-kl_sudu)*abs(qx_sudu-kl_sudu)
     &	/4/kl_zhijing/kl_midu-2*kl_r_sudu*(kl_sudu/ban_jing+jiao_sudu)
	elseif(re<=1000.0) then
		cd=98.33/re-2778/re/re+0.3644
		dvdtc=3*cd*qx_midu*(qx_sudu-kl_sudu)*abs(qx_sudu-kl_sudu)
     &	/4/kl_zhijing/kl_midu-2*kl_r_sudu*(kl_sudu/ban_jing+jiao_sudu)
	elseif(re<=5000.0) then
		cd=148.62/re-47500/re/re+0.357
		dvdtc=3*cd*qx_midu*(qx_sudu-kl_sudu)*abs(qx_sudu-kl_sudu)
     &	/4/kl_zhijing/kl_midu-2*kl_r_sudu*(kl_sudu/ban_jing+jiao_sudu)
	elseif(re<=10000.0) then
		cd=-490.546/re+578700/re/re+0.46
		dvdtc=3*cd*qx_midu*(qx_sudu-kl_sudu)*abs(qx_sudu-kl_sudu)
     &	/4/kl_zhijing/kl_midu-2*kl_r_sudu*(kl_sudu/ban_jing+jiao_sudu)
	elseif(re<=50000.0) then
		cd=-1662.5/re+5416700/re/re+0.5191
		dvdtc=3*cd*qx_midu*(qx_sudu-kl_sudu)*abs(qx_sudu-kl_sudu)
     &	/4/kl_zhijing/kl_midu-2*kl_r_sudu*(kl_sudu/ban_jing+jiao_sudu)
	endif
	end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	����Һ�����򵱵�ʱ�䲽��
	function delta_c(a,b,c,d,e,f,cvr,cr,co)
	double precision a,b,c,d,e,f,cvr,cr,co
	double precision qx_midu,qx_sudu,kl_midu,kl_sudu,kl_zhijing,
     &	qx_nianxing,kl_r_sudu,ban_jing,jiao_sudu
	double precision delta_c
	qx_midu=a
	qx_sudu=b
	kl_midu=c
	kl_sudu=d
	kl_zhijing=e
	qx_nianxing=f
	kl_r_sudu=cvr
	ban_jing=cr
	jiao_sudu=co
	if (abs(qx_sudu-kl_sudu)<1D-8) then
c		delta_c=(qx_sudu-kl_sudu)/(18*qx_nianxing*(qx_sudu-kl_sudu)/
c     &	kl_zhijing/kl_zhijing/kl_midu+2*kl_r_sudu*(kl_sudu/ban_jing+
c     &	jiao_sudu))
c		delta_c=1.0/(18*qx_nianxing/kl_zhijing/kl_zhijing/kl_midu)
		delta_c=kl_zhijing*kl_zhijing*kl_midu/18/qx_nianxing
	else
c		delta_c=1.0/dvdtc(qx_midu,qx_sudu,kl_midu,kl_sudu,kl_zhijing,
c     &	qx_nianxing,kl_r_sudu,ban_jing,jiao_sudu)*(qx_sudu-kl_sudu)
c		delta_c=(qx_sudu-kl_sudu)/dvdtc(qx_midu,qx_sudu,kl_midu,
c     &	kl_sudu,kl_zhijing,qx_nianxing,kl_r_sudu,ban_jing,jiao_sudu)
		delta_c=(qx_sudu-kl_sudu)/dvdtc(qx_midu,qx_sudu,kl_midu,
     &	kl_sudu,kl_zhijing,qx_nianxing,0D0,ban_jing,jiao_sudu)
	endif
	delta_c=delta_c*0.3
	end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	����Һ�ξ����ٶ�ƫ����
	function dvdtr(a,b,c,d,e,f,rvc,rr,ro)
	double precision a,b,c,d,e,f,rvc,rr,ro
	double precision re,cd,qx_midu,qx_sudu,kl_midu,kl_sudu,kl_zhijing,
     &	qx_nianxing,kl_c_sudu,ban_jing,jiao_sudu,dvdtr
	qx_midu=a
	qx_sudu=b
	kl_midu=c
	kl_sudu=d
	kl_zhijing=e
	qx_nianxing=f
	kl_c_sudu=rvc
	ban_jing=rr
	jiao_sudu=ro
	re=qx_midu*abs(qx_sudu-kl_sudu)*kl_zhijing/qx_nianxing
	if (re<1D-8) then
		dvdtr=18*qx_nianxing*(qx_sudu-kl_sudu)/kl_zhijing/kl_zhijing/
     &	kl_midu+ban_jing*(kl_c_sudu/ban_jing+jiao_sudu)**2
	elseif (re<0.1) then
		cd=24.0/re
		dvdtr=3*cd*qx_midu*(qx_sudu-kl_sudu)*abs(qx_sudu-kl_sudu)/4/
     &	kl_zhijing/kl_midu+ban_jing*(kl_c_sudu/ban_jing+jiao_sudu)**2
	elseif(re<=1.0) then
		cd=22.73/re+0.0903/re/re+3.69
		dvdtr=3*cd*qx_midu*(qx_sudu-kl_sudu)*abs(qx_sudu-kl_sudu)/4/
     &	kl_zhijing/kl_midu+ban_jing*(kl_c_sudu/ban_jing+jiao_sudu)**2
	elseif(re<=10.0) then
		cd=29.1667/re-3.8889/re/re+1.222
		dvdtr=3*cd*qx_midu*(qx_sudu-kl_sudu)*abs(qx_sudu-kl_sudu)/4/
     &	kl_zhijing/kl_midu+ban_jing*(kl_c_sudu/ban_jing+jiao_sudu)**2
	elseif(re<=100.0) then
		cd=46.5/re-116.67/re/re+0.6167
		dvdtr=3*cd*qx_midu*(qx_sudu-kl_sudu)*abs(qx_sudu-kl_sudu)/4/
     &	kl_zhijing/kl_midu+ban_jing*(kl_c_sudu/ban_jing+jiao_sudu)**2
	elseif(re<=1000.0) then
		cd=98.33/re-2778/re/re+0.3644
		dvdtr=3*cd*qx_midu*(qx_sudu-kl_sudu)*abs(qx_sudu-kl_sudu)/4/
     &	kl_zhijing/kl_midu+ban_jing*(kl_c_sudu/ban_jing+jiao_sudu)**2
	elseif(re<=5000.0) then
		cd=148.62/re-47500/re/re+0.357
		dvdtr=3*cd*qx_midu*(qx_sudu-kl_sudu)*abs(qx_sudu-kl_sudu)/4/
     &	kl_zhijing/kl_midu+ban_jing*(kl_c_sudu/ban_jing+jiao_sudu)**2
	elseif(re<=10000.0) then
		cd=-490.546/re+578700/re/re+0.46
		dvdtr=3*cd*qx_midu*(qx_sudu-kl_sudu)*abs(qx_sudu-kl_sudu)/4/
     &	kl_zhijing/kl_midu+ban_jing*(kl_c_sudu/ban_jing+jiao_sudu)**2
	elseif(re<=50000.0) then
		cd=-1662.5/re+5416700/re/re+0.5191
		dvdtr=3*cd*qx_midu*(qx_sudu-kl_sudu)*abs(qx_sudu-kl_sudu)/4/
     &	kl_zhijing/kl_midu+ban_jing*(kl_c_sudu/ban_jing+jiao_sudu)**2
	endif
	end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	����Һ�ξ��򵱵�ʱ�䲽��
	function delta_r(a,b,c,d,e,f,rvc,rr,ro)
	double precision a,b,c,d,e,f,rvc,rr,ro
	double precision qx_midu,qx_sudu,kl_midu,kl_sudu,kl_zhijing,
     &	qx_nianxing,kl_c_sudu,ban_jing,jiao_sudu
	double precision delta_r
	qx_midu=a
	qx_sudu=b
	kl_midu=c
	kl_sudu=d
	kl_zhijing=e
	qx_nianxing=f
	kl_c_sudu=rvc
	ban_jing=rr
	jiao_sudu=ro
	if (abs(qx_sudu-kl_sudu)<1D-10) then
c		delta_r=1.0/(18*qx_nianxing/kl_zhijing/kl_zhijing/kl_midu+
c     &	ban_jing*(kl_c_sudu/ban_jing+jiao_sudu)**2/(qx_sudu-kl_sudu))
		delta_r=kl_zhijing*kl_zhijing*kl_midu/18/qx_nianxing
	else
c		delta_r=1.0/dvdtr(qx_midu,qx_sudu,kl_midu,kl_sudu,kl_zhijing,
c     &	qx_nianxing,kl_c_sudu,ban_jing,jiao_sudu)*(qx_sudu-kl_sudu)
		delta_r=(qx_sudu-kl_sudu)/dvdtr(qx_midu,qx_sudu,kl_midu,
     &	kl_sudu,kl_zhijing,qx_nianxing,0D0,ban_jing,0D0)
	endif
	delta_r=delta_r*0.5
c	delta_r=delta_r*2.7
	end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	�����������������㺬������һԪ����
	subroutine fsim(a,b,ns,one_function,g)
	external one_function
	double precision a,b,one_function,g
	double precision at,bt,ht,st,fat,fbt,xt,f2t,f4t
	integer ns,nt
	at=a
	bt=b
	nt=20
	ht=(bt-at)/nt/2.0
	st=0.0
	fat=one_function(at,ns)
	fbt=one_function(bt,ns)
	xt=at+ht
	f2t=0.0
	f4t=one_function(xt,ns)
	do j=1,nt-1
		xt=xt+ht
		f2t=f2t+one_function(xt,ns)
		xt=xt+ht
		f4t=f4t+one_function(xt,ns)
	enddo
	g=ht/3.0*(fat+fbt+4.0*f4t+2.0*f2t)
	end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	��������ѧ��ʧ���ֺ���
	function thermo_loss(x,ns)
	use WASPCN
	implicit double precision (a-h,o-z)
	parameter (NStage=10,NSpan=11,NStream=11)
	common/com13/static_pressure(0:2*NStage*NStream-1)
	common/com14/static_temperature(0:2*NStage*NStream-1)
	common/com15/d_mass_fraction(0:2*NStage*NStream-1)
	common/com16/axial_distance(0:2*NStage*NStream-1)
	common/com17/a_mass_flowrate(0:2*NStage*NStream-1)
	double precision wetness_massflow(0:2*NStage*NStream-1)
	double precision sp,st,sat
	double precision hh_v,hh_l,tt_s
	double precision latent_heat,dqdx,super_cooling
	double precision temp
	integer ns,ntemp
!	sp��ѹ��st���£�sat�����¶�
!	hh_v����������ֵ��hh_l����ˮ��ֵ��tt_s�����¶�
!	�����Բ�ֵ��������������ֵ
	do i=0,2*NStage*NStream-2
	if (axial_distance(i)>axial_distance(i+1)) then
		temp=axial_distance(i)
		axial_distance(i)=axial_distance(i+1)
		axial_distance(i+1)=temp
!
		temp=static_pressure(i)
		static_pressure(i)=static_pressure(i+1)
		static_pressure(i+1)=temp
!
		temp=static_temperature(i)
		static_temperature(i)=static_temperature(i+1)
		static_temperature(i+1)=temp
!
		temp=d_mass_fraction(i)
		d_mass_fraction(i)=d_mass_fraction(i+1)
		d_mass_fraction(i+1)=temp
!
		temp=a_mass_flowrate(i)
		a_mass_flowrate(i)=a_mass_flowrate(i+1)
		a_mass_flowrate(i+1)=temp
	endif
	enddo
!	������������ַǵ������������������λ��
!	Ǳ�ڵ����⣬����Ҫ��������Դ!
	do i=0,2*NStage*NStream-1
		wetness_massflow(i)=d_mass_fraction(i)*a_mass_flowrate(i)
	enddo
	call SPLINE(axial_distance,wetness_massflow,2*NStage*NStream-1,
     &	x,temp,dqdx)
	do i=0,2*NStage*NStream-2
	if (x>=axial_distance(i).and.x<=axial_distance(i+1)) then
	sp=static_pressure(i)+(static_pressure(i+1)-static_pressure(i))/
     &	(axial_distance(i+1)-axial_distance(i))*(x-axial_distance(i))
	st=static_temperature(i)+
     &	(static_temperature(i+1)-static_temperature(i))/
     &	(axial_distance(i+1)-axial_distance(i))*(x-axial_distance(i))
	endif
	enddo
	call P2HL(sp/1000000,hh_l,ntemp)
	call P2HG(sp/1000000,hh_v,ntemp)
	latent_heat=(hh_v-hh_l)*1000
	call P2T(sp/1000000,tt_s,ntemp)
	tt_s=tt_s+273.15
	if (tt_s>=st) then
		super_cooling=tt_s-st
	else
		super_cooling=0.0
	endif
	sat=tt_s
	if (dqdx<=0.0) dqdx=0.0
	thermo_loss=latent_heat*dqdx*super_cooling/sat
	end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	��������������������ػ���
	SUBROUTINE FSIM2(A,B,D_L_BOUNDARY,DRAG_LOSS,L,G)
	EXTERNAL DRAG_LOSS,D_L_BOUNDARY
	double precision A,B,DRAG_LOSS,G
	double precision H,X,S,FA,FB,F2,F4,G1,G2
	INTEGER N,L
	N=50
	H=(B-A)/N/2.0
	S=0.0
	CALL FSIM1(A,D_L_BOUNDARY,DRAG_LOSS,L,FA)
	CALL FSIM1(B,D_L_BOUNDARY,DRAG_LOSS,L,FB)
	X=A+H
	F2=0.0
	CALL FSIM1(X,D_L_BOUNDARY,DRAG_LOSS,L,F4)
	DO J=1,N-1
		X=X+H
	CALL FSIM1(X,D_L_BOUNDARY,DRAG_LOSS,L,G1)
		F2=F2+G1
		X=X+H
	CALL FSIM1(X,D_L_BOUNDARY,DRAG_LOSS,L,G2)
		F4=F4+G2
	ENDDO
	G=H/3.0*(FA+FB+4.0*F4+2.0*F2)
	END
!	������������������һԪ����
	SUBROUTINE FSIM1(X,D_L_BOUNDARY,DRAG_LOSS,L,S)
	EXTERNAL DRAG_LOSS,D_L_BOUNDARY
	double precision X,DRAG_LOSS,S
	double precision A,B,H,XT,F2,F4,FA,FB
	INTEGER N,L
	CALL D_L_BOUNDARY(X,A,B)
	N=20
	H=(B-A)/N/2.0
	S=0.0
	FA=DRAG_LOSS(X,A,L)
	FB=DRAG_LOSS(X,B,L)
	XT=A+H
	F2=0.0
	F4=DRAG_LOSS(X,XT,L)
	DO J=1,N-1
		XT=XT+H
		F2=F2+DRAG_LOSS(X,XT,L)
		XT=XT+H
		F4=F4+DRAG_LOSS(X,XT,L)
	ENDDO
	S=H/3.0*(FA+FB+4.0*F4+2.0*F2)
	END
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	����droplet re-entrainment��ʧ���ֱ��ʽ
	FUNCTION DRAG_LOSS(x0,y0,L)
	implicit double precision (a-h,o-z)
	parameter (NStage=10,NSpan=11,NStream=11)
	common/com01/cylindrical_r(0:2*NStage*NStream-1,0:NSpan-1)
	common/com02/cylindrical_z(0:2*NStage*NStream-1,0:NSpan-1)
	common/com04/velocity_axial(0:2*NStage*NStream-1,0:NSpan-1)
	common/com05/velocity_radial(0:2*NStage*NStream-1,0:NSpan-1)
	common/com06/velocity_circumferential(0:2*NStage*NStream-1,
     &	0:NSpan-1)
	common/com07/density_vapor(0:2*NStage*NStream-1,0:NSpan-1)
	common/com09/particle_diameter(0:2*NStage*NStream-1,0:NSpan-1)
	common/com10/dynamic_viscosity(0:2*NStage*NStream-1,0:NSpan-1)
	common/com18/volume_fraction(0:2*NStage*NStream-1,0:NSpan-1)
	common/com30/v_d_a(0:2*NStage*NStream-1,0:NSpan-1)
	common/com31/v_d_r(0:2*NStage*NStream-1,0:NSpan-1)
	common/com32/v_d_c(0:2*NStage*NStream-1,0:NSpan-1)
	double precision DRAG_A_LOSS,x0,y0
!	x0:����,y0:����
	INTEGER L
!	L=1Ϊ����L=2Ϊ����L=3Ϊ����
	double precision VV,VD,RHOV,PDD,DVV,VFD
	double precision re,cd
	double precision x(0:2*NStage*NStream-1,0:NSpan-1)
	double precision y(0:2*NStage*NStream-1,0:NSpan-1)
!	VV�����ٶȣ�VDҺ���ٶȣ�RHOV�����ܶȣ�PDDҺ��ֱ����DVV����ճ�ԣ�VFDҺ���������
!	�����Բ�ֵ��������������ֵ
	do i=0,2*NStage*NStream-1
		do j=0,NSpan-1
			x(i,j)=cylindrical_z(i,j)
			y(i,j)=cylindrical_r(i,j)
		enddo
	enddo
	do 10 i=0,2*NStage*NStream-2
		do 20 j=0,NSpan-2
	if ((((y(i+1,j)-y(i,j))*x0+(x(i,j)-x(i+1,j))*y0+
     &	x(i+1,j)*y(i,j)-x(i,j)*y(i+1,j))*
     &	 ((y(i+1,j+1)-y(i,j+1))*x0+(x(i,j+1)-x(i+1,j+1))*y0+
     &	x(i+1,j+1)*y(i,j+1)-x(i,j+1)*y(i+1,j+1))
     &	<=0.0) .and. 
     &	(((y(i,j+1)-y(i,j))*x0+(x(i,j)-x(i,j+1))*y0+
     &	x(i,j+1)*y(i,j)-x(i,j)*y(i,j+1))*
     &	((y(i+1,j+1)-y(i+1,j))*x0+(x(i+1,j)-x(i+1,j+1))*y0+
     &	x(i+1,j+1)*y(i+1,j)-x(i+1,j)*y(i+1,j+1))
     &	<=0.0)) then
		if (L==1) then
	VV=(velocity_axial(i,j)*(x(i+1,j)-x0)*(y(i,j+1)-y0)+
     &    velocity_axial(i+1,j)*(x0-x(i,j))*(y(i,j+1)-y0)+
     &    velocity_axial(i,j+1)*(x(i+1,j)-x0)*(y0-y(i,j))+
     &    velocity_axial(i+1,j+1)*(x0-x(i,j))*(y0-y(i,j)))/
     &	(x(i+1,j)-x(i,j))/(y(i,j+1)-y(i,j))
	VD=(v_d_a(i,j)*(x(i+1,j)-x0)*(y(i,j+1)-y0)+
     &    v_d_a(i+1,j)*(x0-x(i,j))*(y(i,j+1)-y0)+
     &    v_d_a(i,j+1)*(x(i+1,j)-x0)*(y0-y(i,j))+
     &    v_d_a(i+1,j+1)*(x0-x(i,j))*(y0-y(i,j)))/
     &	(x(i+1,j)-x(i,j))/(y(i,j+1)-y(i,j))
		elseif (L==2) then
	VV=(velocity_circumferential(i,j)*(x(i+1,j)-x0)*(y(i,j+1)-y0)+
     &    velocity_circumferential(i+1,j)*(x0-x(i,j))*(y(i,j+1)-y0)+
     &    velocity_circumferential(i,j+1)*(x(i+1,j)-x0)*(y0-y(i,j))+
     &    velocity_circumferential(i+1,j+1)*(x0-x(i,j))*(y0-y(i,j)))/
     &	(x(i+1,j)-x(i,j))/(y(i,j+1)-y(i,j))
	VD=(v_d_c(i,j)*(x(i+1,j)-x0)*(y(i,j+1)-y0)+
     &    v_d_c(i+1,j)*(x0-x(i,j))*(y(i,j+1)-y0)+
     &    v_d_c(i,j+1)*(x(i+1,j)-x0)*(y0-y(i,j))+
     &    v_d_c(i+1,j+1)*(x0-x(i,j))*(y0-y(i,j)))/
     &	(x(i+1,j)-x(i,j))/(y(i,j+1)-y(i,j))
		elseif (L==3) then
	VV=(velocity_radial(i,j)*(x(i+1,j)-x0)*(y(i,j+1)-y0)+
     &    velocity_radial(i+1,j)*(x0-x(i,j))*(y(i,j+1)-y0)+
     &    velocity_radial(i,j+1)*(x(i+1,j)-x0)*(y0-y(i,j))+
     &    velocity_radial(i+1,j+1)*(x0-x(i,j))*(y0-y(i,j)))/
     &	(x(i+1,j)-x(i,j))/(y(i,j+1)-y(i,j))
	VD=(v_d_r(i,j)*(x(i+1,j)-x0)*(y(i,j+1)-y0)+
     &    v_d_r(i+1,j)*(x0-x(i,j))*(y(i,j+1)-y0)+
     &    v_d_r(i,j+1)*(x(i+1,j)-x0)*(y0-y(i,j))+
     &    v_d_r(i+1,j+1)*(x0-x(i,j))*(y0-y(i,j)))/
     &	(x(i+1,j)-x(i,j))/(y(i,j+1)-y(i,j))
		endif
	endif
20		continue
10	continue
	call SPLINE_AREA(density_vapor,x0,y0,RHOV)
	call SPLINE_AREA(particle_diameter,x0,y0,PDD)
	call SPLINE_AREA(dynamic_viscosity,x0,y0,DVV)
	call SPLINE_AREA(volume_fraction,x0,y0,VFD)
	re=RHOV*abs(VV-VD)*PDD/DVV
	if (re<1D-8) then
		cd=0.0
	elseif (re<0.1) then
		cd=24.0/re
	elseif(re<=1.0) then
		cd=22.73/re+0.0903/re/re+3.69
	elseif(re<=10.0) then
		cd=29.1667/re-3.8889/re/re+1.222
	elseif(re<=100.0) then
		cd=46.5/re-116.67/re/re+0.6167
	elseif(re<=1000.0) then
		cd=98.33/re-2778/re/re+0.3644
	elseif(re<=5000.0) then
		cd=148.62/re-47500/re/re+0.357
	elseif(re<=10000.0) then
		cd=-490.546/re+578700/re/re+0.46
	elseif(re<=50000.0) then
		cd=-1662.5/re+5416700/re/re+0.5191
	endif
	DRAG_LOSS=(3.0*3.1415926/2.0/PDD)*cd*RHOV*abs((VV-VD)**3)*VFD*y0
c	DRAG_LOSS=(3.0*3.1415926/2.0/PDD)*cd*RHOV*((VV-VD)**3)*VFD*y0
	END
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	����droplet re-entrainment��ʧ�����ޱ��ʽ
	SUBROUTINE D_L_BOUNDARY(x,y1,y2)
	implicit double precision (a-h,o-z)
	parameter (NStage=10,NSpan=11,NStream=11)
	common/com01/cylindrical_r(0:2*NStage*NStream-1,0:NSpan-1)
	common/com02/cylindrical_z(0:2*NStage*NStream-1,0:NSpan-1)
	double precision x,y1,y2
!	y1���ޣ�y2����
	do i=0,2*NStage*NStream-2
	if (x>=cylindrical_z(i,0).and.
     &	x<=cylindrical_z(i+1,0)) then
	y1=cylindrical_r(i,0)+(cylindrical_r(i+1,0)-cylindrical_r(i,0))/
     &  (cylindrical_z(i+1,0)-cylindrical_z(i,0))*(x-cylindrical_z(i,0))
	endif
	if (x>=cylindrical_z(i,NSpan-1).and.
     &	x<=cylindrical_z(i+1,NSpan-1)) then
	y2=cylindrical_r(i,NSpan-1)+(cylindrical_r(i+1,NSpan-1)-
     &	cylindrical_r(i,NSpan-1))/(cylindrical_z(i+1,NSpan-1)-
     &	cylindrical_z(i,NSpan-1))*(x-cylindrical_z(i,NSpan-1))
	endif
	enddo
	END
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	����Droplet impact��ʧ���ֱ��ʽ
!	��������ȡΪ��ҶǰԵ����ֵ
	FUNCTION DEOPLET_IMPACT(r,NS)
	implicit double precision (a-h,o-z)
	parameter (NStage=10,NSpan=11,NStream=11)
	parameter (zhuansu=1500*2*3.1415926/60.0)
	common/com02/cylindrical_z(0:2*NStage*NStream-1,0:NSpan-1)
	common/com08/density_liquid(0:2*NStage*NStream-1,0:NSpan-1)
	common/com18/volume_fraction(0:2*NStage*NStream-1,0:NSpan-1)
	common/com19/sound_speed(0:2*NStage*NStream-1,0:NSpan-1)
	common/com33/v_ds_a(0:2*NStage*NStream-1,0:NSpan-1)
	common/com34/v_ds_r(0:2*NStage*NStream-1,0:NSpan-1)
	common/com35/v_ds_c(0:2*NStage*NStream-1,0:NSpan-1)
	common/com41/Number_r(NStage)
	common/com56/idirection
	double precision r
!	rȡΪ��ҶǰԵ�ľ����������
	integer NS
	double precision droplet_density,droplet_sound,droplet_fraction
	double precision c1a,c1c
!	c1aΪ����Һ�������ٶ�,c1cΪ����Һ�������ٶ�
	double precision lr,pit
	double precision steel_density,steel_sound
	double precision w,p,temp
	double precision tijifenshu(0:NSpan-1)
	double precision slr(NStage,0:NSpan-1),slz(NStage,0:NSpan-1)
	double precision str(NStage,0:NSpan-1),stz(NStage,0:NSpan-1)
	double precision rlr(NStage,0:NSpan-1),rlz(NStage,0:NSpan-1)
	double precision rtr(NStage,0:NSpan-1),rtz(NStage,0:NSpan-1)
	double precision kl_s_weizhi(0:NStream-1)
	double precision kl_midu(0:NStream-1)
	double precision kl_shengsu(0:NStream-1)
	double precision kl_tijifenshu(0:NStream-1)
	double precision kl_r_weizhi(0:NStream-1)
	double precision kl_asudu(0:NStream-1)
	double precision kl_csudu(0:NStream-1)
	double precision kldmidu(NStage,0:NSpan-1)
	double precision kldshengsu(NStage,0:NSpan-1)
	double precision kldtijifenshu(NStage,0:NSpan-1)
	double precision kldasudu(NStage,0:NSpan-1)
	double precision kldcsudu(NStage,0:NSpan-1)
	double precision rr(0:NSpan-1)
	double precision midu(0:NSpan-1)
	double precision shengsu(0:NSpan-1)
	double precision asudu(0:NSpan-1)
	double precision csudu(0:NSpan-1)
	double precision uu(0:NSpan-1)
	double precision wsudu(0:NSpan-1)
	double precision beta(0:NSpan-1)
	double precision length(0:NSpan-1)
	double precision pitch(0:NSpan-1)
	steel_density=7.85*1000.0
	steel_sound=5920.0
	call stator_rotor_a_r(slr,slz,str,stz,rlr,rlz,rtr,rtz)
	do k=0,NSpan-1
	do j=0,NStream-1
		kl_s_weizhi(j)=cylindrical_z(2*(NS-1)*NStream+j,k)
		kl_midu(j)=density_liquid(2*(NS-1)*NStream+j,k)
		kl_shengsu(j)=sound_speed(2*(NS-1)*NStream+j,k)
		kl_tijifenshu(j)=volume_fraction(2*(NS-1)*NStream+j,k)
!	��Ҷ�������
		kl_r_weizhi(j)=cylindrical_z((2*NS-1)*NStream+j,k)
		kl_asudu(j)=v_ds_a((2*NS-1)*NStream+j,k)
		kl_csudu(j)=v_ds_c((2*NS-1)*NStream+j,k)
!	��Ҷ�������
	enddo
	call SPLINE(kl_s_weizhi,kl_midu,NStream-1,stz(NS,k),
     &	kldmidu(NS,k),temp)
	call SPLINE(kl_s_weizhi,kl_shengsu,NStream-1,stz(NS,k),
     &	kldshengsu(NS,k),temp)
	call SPLINE(kl_s_weizhi,kl_tijifenshu,NStream-1,stz(NS,k),
     &	kldtijifenshu(NS,k),temp)
!	��ҶβԵ����
	call SPLINE(kl_r_weizhi,kl_asudu,NStream-1,rlz(NS,k),
     &	kldasudu(NS,k),temp)
	call SPLINE(kl_r_weizhi,kl_csudu,NStream-1,rlz(NS,k),
     &	kldcsudu(NS,k),temp)
!	��ҶǰԵ����
	enddo
	do i=0,NSpan-1
		rr(i)=rlr(NS,i)
!	��ҶǰԵҶ��
		midu(i)=kldmidu(NS,i)
		shengsu(i)=kldshengsu(NS,i)
		tijifenshu(i)=kldtijifenshu(NS,i)
!	��ҶβԵ����		
		asudu(i)=kldasudu(NS,i)
		csudu(i)=kldcsudu(NS,i)
		uu(i)=zhuansu*rlr(NS,i)
		wsudu(i)=sqrt(asudu(i)**2+(uu(i)-csudu(i))**2)
		beta(i)=asin(asudu(i)/wsudu(i))
!	��ҶǰԵ����
	enddo
	call impact_zone(NS,beta,length,pitch)
	call SPLINE(rr,midu,NSpan-1,r,droplet_density,temp)
	call SPLINE(rr,shengsu,NSpan-1,r,droplet_sound,temp)
	call SPLINE(rr,tijifenshu,NSpan-1,r,droplet_fraction,temp)
!	��ҶβԵ����
	call SPLINE(rr,wsudu,NSpan-1,r,w,temp)
	call SPLINE(rr,length,NSpan-1,r,lr,temp)
	call SPLINE(rr,pitch,NSpan-1,r,pit,temp)
!	����Һ�������������������
!	�ο���ʪ�������о�����ͺ˵����ֻ��ڱ���ʪ�����о�		
!	����ƶ��Ƕ�Ϊ�۽ǣ�����ײ������Ҷѹ���棬��Ҫ�ж�
	p=droplet_density*droplet_sound*w/(1+(droplet_density*
     &	droplet_sound)/(steel_density*steel_sound))
	DEOPLET_IMPACT=Number_r(NS)*zhuansu*p*lr*pit*r*droplet_fraction
	END
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	�����ƶ������դ��ٷ���
	SUBROUTINE impact_zone(NS,beta,length,pitch)
	implicit double precision (a-h,o-z)
	parameter (NStage=10,NSpan=11,NStream=11)
	parameter (zhuansu=1500*2*3.1415926/60.0)
	parameter (MS=100,MN=800)
	common/com41/Number_r(NStage)
	common/com56/idirection
	integer NS
	double precision beta(0:NSpan-1)
!	beta:�ƶ��Ƕ�
	double precision beta1(0:NSpan-1),beta2(0:NSpan-1)
!	beta1:�ƶ��Ƕ����,beta2:�ƶ��Ƿ����-beta1
	double precision x0,y0,xp,yp,xn,yn,xm,ym
!	��ת��������,ȡΪ��ҶǰԵ��
!	��Ҷ��������,ȡΪ��ҶǰԵ��
!	��Ҷ������ת����,ȡΪ��ҶǰԵ��
!	��������
	double precision length(0:NSpan-1)
	double precision pitch(0:NSpan-1)
!	length:�����ƶ�����,pitch:դ��ٷֱ�
	integer NSNumber_rns(NStage),r_number
	integer istart,iend
	double precision r(NStage,0:NSpan-1,0:MN)
	double precision z(NStage,0:NSpan-1,0:MN)
	double precision c(NStage,0:NSpan-1,0:MN)
	double precision x1(0:MN),y1(0:MN)
	double precision x2(0:MN),y2(0:MN)
	double precision x3(0:MN),y3(0:MN)
	double precision xtemp,ytemp,temp
!!!!!!
c	print *,'���Ǽ��붯Ҷ��ת�������idirection'
	call cylinder_r_rcz(NSNumber_rns,z,r,c)
	do i=0,NSpan-1
	x0=c(NS,i,0)
	y0=z(NS,i,0)
	xp=c(NS,i,0)+r(NS,i,0)*(2.0*3.1415926/Number_r(NS))
	yp=z(NS,i,0)
!!!!!!
	beta1(i)=3.1415926/2-beta(i)
	beta2(i)=-beta1(i)
		do j=0,NSNumber_rns(NS)-1
			x1(j)=c(NS,i,j)
			y1(j)=z(NS,i,j)
		enddo
		do j=0,NSNumber_rns(NS)-1
			x2(j)=(x1(j)-x0)*cos(beta2(i))-(y1(j)-y0)*sin(beta2(i))+x0
			y2(j)=(x1(j)-x0)*sin(beta2(i))+(y1(j)-y0)*cos(beta2(i))+y0
		enddo
		xn=(xp-x0)*cos(beta2(i))-(yp-y0)*sin(beta2(i))+x0
		yn=(xp-x0)*sin(beta2(i))+(yp-y0)*cos(beta2(i))+y0
!	������ת��Ҷ������������
		do j=0,NSNumber_rns(NS)-2
			if (x2(j)<x2(j+1)) then
				iend=j+1
			endif
		enddo
		do j=NSNumber_rns(NS)-2,0,-1
			if (x2(j)<x2(j+1)) then
				istart=j		
			endif 				
		enddo
		r_number=iend-istart+1
		do j=0,r_number-1
			x3(j)=x2(istart+j)
			y3(j)=y2(istart+j)
		enddo
!	������ת�󵥵�������
		if (xn<=x3(r_number-1)) then
			xtemp=xn
			call SPLINE(x3,y3,r_number-1,xtemp,ytemp,temp)
			xm=(xtemp-x0)*cos(beta1(i))-(ytemp-y0)*sin(beta1(i))+x0
			ym=(xtemp-x0)*sin(beta1(i))+(ytemp-y0)*cos(beta1(i))+y0
			length(i)=ym-y0
			pitch(i)=1.0
		else
			xtemp=x3(r_number-1)
			ytemp=y3(r_number-1)
			xm=(xtemp-x0)*cos(beta1(i))-(ytemp-y0)*sin(beta1(i))+x0
			ym=(xtemp-x0)*sin(beta1(i))+(ytemp-y0)*cos(beta1(i))+y0
			length(i)=ym-y0
	 	pitch(i)=(xm-x0)/(r(NS,i,0)*(2.0*3.1415926/Number_r(NS)))
c	�ƶ���ƫ��ʱ���ھ�ٷֱȻ�������⣬��Ҫ������
		if (pitch(i)<0.0) then
	print *,'the impact angle or the blade profile have some problem!'	
		endif
		endif
!	���㽻������
!	xtemp,ytempΪ��ת��Ľ�����������������ֵ
!	xm,ymΪԭʼ�Ľ���������������е�����
	enddo
	END
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	���㾲Ҷ�Ͷ�Ҷǰ��Ե����;�������
	SUBROUTINE stator_rotor_a_r(slr,slz,str,stz,rlr,rlz,rtr,rtz)
	implicit double precision (a-h,o-z)
	parameter (NStage=10,NSpan=11,NStream=11)
	parameter (MS=100,MN=800)
	common/com01/cylindrical_r(0:2*NStage*NStream-1,0:NSpan-1)
	common/com02/cylindrical_z(0:2*NStage*NStream-1,0:NSpan-1)
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
	common/com60/rss(NStage,0:MS,0:MN)
	common/com61/tss(NStage,0:MS,0:MN)
	common/com62/ass(NStage,0:MS,0:MN)
	common/com63/rps(NStage,0:MS,0:MN)
	common/com64/tps(NStage,0:MS,0:MN)
	common/com65/aps(NStage,0:MS,0:MN)
!	��Ҷ�������ѹ����Բ��������Ҷչ������Ƶ�����
	common/com72/rsrn(NStage,0:MS,0:MN)
	common/com73/tsrn(NStage,0:MS,0:MN)
	common/com74/asrn(NStage,0:MS,0:MN)
!	��Ҷ������Բ��������Ҷչ������Ƶ��������
	double precision stator_lr(NStage,0:MS),stator_lz(NStage,0:MS)
	double precision stator_tr(NStage,0:MS),stator_tz(NStage,0:MS)
	double precision rotor_lr(NStage,0:MS),rotor_lz(NStage,0:MS)
	double precision rotor_tr(NStage,0:MS),rotor_tz(NStage,0:MS)
!	��Ҷ�Ͷ�Ҷǰ��Ե����������;�������ֵ
	double precision xlz(0:MS),ylr(0:MS)
	double precision xtz(0:MS),ytr(0:MS)
	double precision x1,y1,x2,y2,xl,yl,xr,yr
	double precision slr(NStage,0:NSpan-1),slz(NStage,0:NSpan-1)
	double precision str(NStage,0:NSpan-1),stz(NStage,0:NSpan-1)
	double precision rlr(NStage,0:NSpan-1),rlz(NStage,0:NSpan-1)
	double precision rtr(NStage,0:NSpan-1),rtz(NStage,0:NSpan-1)
!	�����������Ͼ�Ҷ�Ͷ�Ҷǰ��Ե����,������ֵ����
	do k=1,NStage
		do i=0,NSection_s(k)-1
!	��i����������
			stator_lz(k,i)=ass(k,i,0)
			stator_tz(k,i)=ass(k,i,0)
			do j=0,NSNumber_s(k,i)-1
!	��Ҷ����������
				if (ass(k,i,j)<=stator_lz(k,i)) then
					stator_lz(k,i)=ass(k,i,j)
					stator_lr(k,i)=rss(k,i,j)
				endif
				if (ass(k,i,j)>=stator_tz(k,i)) then
					stator_tz(k,i)=ass(k,i,j)
					stator_tr(k,i)=rss(k,i,j)
				endif
			enddo
			do j=0,NPNumber_s(k,i)-1
!	��Ҷѹ��������
				if (aps(k,i,j)<=stator_lz(k,i)) then
					stator_lz(k,i)=aps(k,i,j)
					stator_lr(k,i)=rps(k,i,j)
				endif
				if (aps(k,i,j)>=stator_tz(k,i)) then
					stator_tz(k,i)=aps(k,i,j)
					stator_tr(k,i)=rps(k,i,j)
				endif
			enddo
		enddo
		do i=0,NSection_r(k)-1
!	��i����������
!	ֱ�ӽ���Ҷ�����漫ֵ��Ϊǰ��Ե����
			rotor_lz(k,i)=asrn(k,i,0)
			rotor_lr(k,i)=rsrn(k,i,0)
			rotor_tz(k,i)=asrn(k,i,NSNumber_rn(k,i)-1)
			rotor_tr(k,i)=rsrn(k,i,NSNumber_rn(k,i)-1)
		enddo
	enddo
!	��������������ֵ����NSpan����ǰ��Ե����
	do k=1,NStage
!	��k��
!	��Ҷǰ��Ե����
		do in=0,NSection_s(k)-1
			xlz(in)=stator_lz(k,in)
			ylr(in)=stator_lr(k,in)
			xtz(in)=stator_tz(k,in)
			ytr(in)=stator_tr(k,in)
		enddo
		do is=0,NSpan-1
			do j=2*(k-1)*NStream,(2*k-1)*NStream-2
				x1=cylindrical_z(j+1,is)
				y1=cylindrical_r(j+1,is)
				x2=cylindrical_z(j,is)
				y2=cylindrical_r(j,is)
				call SPLINE_INTERSECTION(ylr,xlz,NSection_s(k)-1,
     &				y1,x1,y2,x2,yl,xl)
				if ((xl-x1)*(xl-x2)<=0.0) then
!	��ҶǰԵ����
					slz(k,is)=xl
					slr(k,is)=yl
				endif
				call SPLINE_INTERSECTION(ytr,xtz,NSection_s(k)-1,
     &				y1,x1,y2,x2,yr,xr)
				if ((xr-x1)*(xr-x2)<=0.0) then
!	��ҶβԵ����
					stz(k,is)=xr
					str(k,is)=yr
				endif
			enddo
		enddo
!	��Ҷǰ��Ե����
		do in=0,NSection_r(k)-1
			xlz(in)=rotor_lz(k,in)
			ylr(in)=rotor_lr(k,in)
			xtz(in)=rotor_tz(k,in)
			ytr(in)=rotor_tr(k,in)
		enddo
		do is=0,NSpan-1
			do j=(2*k-1)*NStream,2*k*NStream-2
				x1=cylindrical_z(j+1,is)
				y1=cylindrical_r(j+1,is)
				x2=cylindrical_z(j,is)
				y2=cylindrical_r(j,is)
				call SPLINE_INTERSECTION(ylr,xlz,NSection_r(k)-1,
     &				y1,x1,y2,x2,yl,xl)
				if ((xl-x1)*(xl-x2)<=0.0) then
!	��ҶǰԵ����
					rlz(k,is)=xl
					rlr(k,is)=yl
				endif
				call SPLINE_INTERSECTION(ytr,xtz,NSection_r(k)-1,
     &				y1,x1,y2,x2,yr,xr)
				if ((xr-x1)*(xr-x2)<=0.0) then
!	��ҶβԵ����
					rtz(k,is)=xr
					rtr(k,is)=yr
				endif
			enddo
		enddo
	enddo
	END
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	���붯Ҷ���������򣬾�����������������ض�Ҷ�����򣬾�����������
	SUBROUTINE cylinder_r_rcz(NSNumber_rns,z,r,c)
	implicit double precision (a-h,o-z)
	parameter (NStage=10,NSpan=11,NStream=11)
	parameter (MS=100,MN=800)
	common/com01/cylindrical_r(0:2*NStage*NStream-1,0:NSpan-1)
	common/com02/cylindrical_z(0:2*NStage*NStream-1,0:NSpan-1)
	common/com51/NSection_r(NStage)
	common/com57/NSNumber_rn(NStage,0:MS)
	common/com72/rsrn(NStage,0:MS,0:MN)
	common/com73/tsrn(NStage,0:MS,0:MN)
	common/com74/asrn(NStage,0:MS,0:MN)
	double precision rr(0:MN)
	double precision zz(0:MN)
	double precision cc(0:MN)
	double precision r(NStage,0:NSpan-1,0:MN)
	double precision z(NStage,0:NSpan-1,0:MN)
	double precision c(NStage,0:NSpan-1,0:MN)
	double precision x1,y1,x2,y2,x3,y3
	double precision temp
	integer NSNumber_rns(NStage)
	do k=1,NStage
		NSNumber_rns(k)=NSNumber_rn(k,0)
		do i=0,NSection_r(k)-1
			if (NSNumber_rn(k,i)<=NSNumber_rns(k)) then
				NSNumber_rns(k)=NSNumber_rn(k,i)
			endif
		enddo
	enddo
!	���������Ҷ����������С���Ƶ����
	do k=1,NStage
	do j=0,NSNumber_rns(k)-1
	do i=0,NSection_r(k)-1
		rr(i)=rsrn(k,i,j)
		zz(i)=asrn(k,i,j)
		cc(i)=tsrn(k,i,j)
	enddo
	do jj=0,NSpan-1
!	��������
	do ii=(2*k-1)*NStream,2*k*NStream-2
	x1=cylindrical_z(ii,jj)
	y1=cylindrical_r(ii,jj)
	x2=cylindrical_z(ii+1,jj)
	y2=cylindrical_r(ii+1,jj)
	call SPLINE_INTERSECTION(rr,zz,NSection_r(k)-1,y2,x2,y1,x1,y3,x3)
	if ((x3-x1)*(x3-x2)<=0.0) then
		z(k,jj,j)=x3
		r(k,jj,j)=y3
	endif
	enddo
!	����������
	y1=r(k,jj,j)
	call SPLINE(rr,cc,NSection_r(k)-1,y1,x1,temp)
	c(k,jj,j)=x1
!!!
	enddo
	enddo
	enddo
	END	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	���������Һ���ٶȲ������Լ����������
!	�����Ӻ������������Ҷ�Ͷ�Ҷǰ��Ե���꣬
!	��Ҷդͨ���ڼ��������Һ������ߣ�
!	���������Һ������߲��죬����Ҷդ�ٷֱ�(������)
c	���������ĳ�����ƫС������������Ҫ���¿��Ǽ��㷽��
	SUBROUTINE DEPOSITION_RATE(dp_rt_s,dp_rt_r)
	parameter (NStage=10,NSpan=11,NStream=11)
	implicit double precision (a-h,o-z)
	common/com01/cylindrical_r(0:2*NStage*NStream-1,0:NSpan-1)
	common/com02/cylindrical_z(0:2*NStage*NStream-1,0:NSpan-1)
	common/com03/omiga(0:2*NStage*NStream-1,0:NSpan-1)
	common/com04/velocity_axial(0:2*NStage*NStream-1,0:NSpan-1)
	common/com05/velocity_radial(0:2*NStage*NStream-1,0:NSpan-1)
	common/com06/velocity_circumferential(0:2*NStage*NStream-1,
     &	0:NSpan-1)
 	common/com30/v_d_a(0:2*NStage*NStream-1,0:NSpan-1)
	common/com31/v_d_r(0:2*NStage*NStream-1,0:NSpan-1)
	common/com32/v_d_c(0:2*NStage*NStream-1,0:NSpan-1)
	common/com40/Number_s(NStage)
	common/com41/Number_r(NStage)
	double precision v_v_a(0:2*NStage*NStream-1,0:NSpan-1)
	double precision v_v_r(0:2*NStage*NStream-1,0:NSpan-1)
	double precision v_v_c(0:2*NStage*NStream-1,0:NSpan-1)
!	�����ٶ�
	double precision v_l_a(0:2*NStage*NStream-1,0:NSpan-1)
	double precision v_l_r(0:2*NStage*NStream-1,0:NSpan-1)
	double precision v_l_c(0:2*NStage*NStream-1,0:NSpan-1)
!	һ��Һ��Һ���ٶ�
	double precision slr(NStage,0:NSpan-1),slz(NStage,0:NSpan-1)
	double precision str(NStage,0:NSpan-1),stz(NStage,0:NSpan-1)
	double precision rlr(NStage,0:NSpan-1),rlz(NStage,0:NSpan-1)
	double precision rtr(NStage,0:NSpan-1),rtz(NStage,0:NSpan-1)
!	���������Ͼ�Ҷ�Ͷ�Ҷǰ��Ե����,������ֵ����
	double precision delta
!	����ʱ�䲽��
	double precision svva(0:100000)
	double precision svpa(0:100000)
	double precision svvc(0:100000)
	double precision svpc(0:100000)
	double precision slva(0:100000)
	double precision slpa(0:100000)
	double precision slvc(0:100000)
	double precision slpc(0:100000)
	double precision rvva(0:100000)
	double precision rvpa(0:100000)
	double precision rvvc(0:100000)
	double precision rvpc(0:100000)
	double precision rlva(0:100000)
	double precision rlpa(0:100000)
	double precision rlvc(0:100000)
	double precision rlpc(0:100000)
!	��Ҷ/��Ҷ,����/Һ��,�ٶ�/λ��,����/����
	double precision weizhi_a(0:NStream-1)
!	Ҷդ��������λ��
	double precision qx_a_sudu(0:NStream-1)
	double precision qx_c_sudu(0:NStream-1)
!	��������������ٶ�
	double precision qx_c_weizhi(NStage,0:NSpan-1)
!	�������������ٶ�
	double precision dp_rt_s(NStage,0:NSpan-1)
!	��Ҷһ��Һ�γ�����
	double precision yx_a_sudu(0:NStream-1)
	double precision yx_c_sudu(0:NStream-1)
!	Һ������������ٶ�
	double precision yx_c_weizhi(NStage,0:NSpan-1)
!	Һ�����������ٶ�
	double precision dp_rt_r(NStage,0:NSpan-1)
!	��Ҷһ��Һ�γ�����
	double precision temp
!	��ʱ����
!!!!!!
!	�������Һ���ٶȸ�ֵ���±�����������������ٶ�ת��Ϊ����ٶ�
	do i=0,2*NStage*NStream-1
		do j=0,NSpan-1
			v_v_a(i,j)=velocity_axial(i,j)
			v_v_r(i,j)=velocity_radial(i,j)
			v_v_c(i,j)=velocity_circumferential(i,j)-
     &			cylindrical_r(i,j)*omiga(i,j)
			v_l_a(i,j)=v_d_a(i,j)
			v_l_r(i,j)=v_d_r(i,j)
			v_l_c(i,j)=v_d_c(i,j)-cylindrical_r(i,j)*omiga(i,j)
		enddo
	enddo
!!!!!!
!	�����Ӻ������㾲Ҷ�Ͷ�Ҷǰ��Ե����
	call stator_rotor_a_r(slr,slz,str,stz,rlr,rlz,rtr,rtz)
!	��ÿһ����Ҷ��Ҷ������м���
	delta=1D-5
	do k=1,NStage
		do j=0,NSpan-1
!!!!!!
!	���㾲Ҷ������
			do i=0,NStream-1
				weizhi_a(i)=cylindrical_z(2*(k-1)*NStream+i,j)
				qx_a_sudu(i)=v_v_a(2*(k-1)*NStream+i,j)
				qx_c_sudu(i)=v_v_c(2*(k-1)*NStream+i,j)
				yx_a_sudu(i)=v_l_a(2*(k-1)*NStream+i,j)
				yx_c_sudu(i)=v_l_c(2*(k-1)*NStream+i,j)
			enddo
			svpa(0)=slz(k,j)
			svpc(0)=0.0
			slpa(0)=slz(k,j)
			slpc(0)=0.0
!	��Ҷͨ���������Һ���ʼ����
			call SPLINE(weizhi_a,qx_a_sudu,NStream-1,
     &			svpa(0),svva(0),temp)
			call SPLINE(weizhi_a,qx_c_sudu,NStream-1,
     &			svpa(0),svvc(0),temp)
			i=1
			do 13 while (svpa(i-1)>=slz(k,j).and.
     &			svpa(i-1)<=stz(k,j))
				call SPLINE(weizhi_a,qx_a_sudu,NStream-1,
     &			svpa(i-1),svva(i-1),temp)
				call SPLINE(weizhi_a,qx_c_sudu,NStream-1,
     &			svpa(i-1),svvc(i-1),temp)
				svpa(i)=svpa(i-1)+svva(i-1)*delta
				svpc(i)=svpc(i-1)+svvc(i-1)*delta
			i=i+1
13			continue
			qx_c_weizhi(k,j)=svpc(i-1)
!	����������������λ��
			call SPLINE(weizhi_a,yx_a_sudu,NStream-1,
     &			slpa(0),slva(0),temp)
			call SPLINE(weizhi_a,yx_c_sudu,NStream-1,
     &			slpa(0),slvc(0),temp)
			i=1
			do 15 while (slpa(i-1)>=slz(k,j).and.
     &			slpa(i-1)<=stz(k,j))
				call SPLINE(weizhi_a,yx_a_sudu,NStream-1,
     &			slpa(i-1),slva(i-1),temp)
				call SPLINE(weizhi_a,yx_c_sudu,NStream-1,
     &			slpa(i-1),slvc(i-1),temp)
				slpa(i)=slpa(i-1)+slva(i-1)*delta
				slpc(i)=slpc(i-1)+slvc(i-1)*delta
			i=i+1
15			continue
			yx_c_weizhi(k,j)=slpc(i-1)
!	����Һ����������λ��
			dp_rt_s(k,j)=abs(qx_c_weizhi(k,j)-yx_c_weizhi(k,j))/
     &			(str(k,j)*(2.0*3.1415926/Number_s(k)))
!!!!!!
!	���㶯Ҷ������
			do i=0,NStream-1
				weizhi_a(i)=cylindrical_z((2*k-1)*NStream+i,j)
				qx_a_sudu(i)=v_v_a((2*k-1)*NStream+i,j)
				qx_c_sudu(i)=v_v_c((2*k-1)*NStream+i,j)
				yx_a_sudu(i)=v_l_a((2*k-1)*NStream+i,j)
				yx_c_sudu(i)=v_l_c((2*k-1)*NStream+i,j)
			enddo
			rvpa(0)=rlz(k,j)
			rvpc(0)=0.0
			rlpa(0)=rlz(k,j)
			rlpc(0)=0.0
!	��Ҷͨ���������Һ���ʼ����
			call SPLINE(weizhi_a,qx_a_sudu,NStream-1,
     &			rvpa(0),rvva(0),temp)
			call SPLINE(weizhi_a,qx_c_sudu,NStream-1,
     &			rvpa(0),rvvc(0),temp)
			i=1
			do 17 while (rvpa(i-1)>=rlz(k,j).and.
     &			rvpa(i-1)<=rtz(k,j))
				call SPLINE(weizhi_a,qx_a_sudu,NStream-1,
     &			rvpa(i-1),rvva(i-1),temp)
				call SPLINE(weizhi_a,qx_c_sudu,NStream-1,
     &			rvpa(i-1),rvvc(i-1),temp)
				rvpa(i)=rvpa(i-1)+rvva(i-1)*delta
				rvpc(i)=rvpc(i-1)+rvvc(i-1)*delta
			i=i+1
17			continue
			qx_c_weizhi(k,j)=rvpc(i-1)
!	����������������λ��
			call SPLINE(weizhi_a,yx_a_sudu,NStream-1,
     &			rlpa(0),rlva(0),temp)
			call SPLINE(weizhi_a,yx_c_sudu,NStream-1,
     &			rlpa(0),rlvc(0),temp)
			i=1
			do 19 while (rlpa(i-1)>=rlz(k,j).and.
     &			rlpa(i-1)<=rtz(k,j))
				call SPLINE(weizhi_a,yx_a_sudu,NStream-1,
     &			rlpa(i-1),rlva(i-1),temp)
				call SPLINE(weizhi_a,yx_c_sudu,NStream-1,
     &			rlpa(i-1),rlvc(i-1),temp)
				rlpa(i)=rlpa(i-1)+rlva(i-1)*delta
				rlpc(i)=rlpc(i-1)+rlvc(i-1)*delta
			i=i+1
19			continue
			yx_c_weizhi(k,j)=rlpc(i-1)
!	����Һ����������λ��
			dp_rt_r(k,j)=abs(qx_c_weizhi(k,j)-yx_c_weizhi(k,j))/
     &			(rtr(k,j)*(2.0*3.1415926/Number_r(k)))
!!!!!!
		enddo
	enddo
!!!!!!
	END
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	�������Ҷդ��Ҷ��Һ������������
	SUBROUTINE FLOWRATE(liquid_s_massflow,liquid_r_massflow,
     &	vapor_s_massflow,vapor_r_massflow)
	implicit double precision (a-h,o-z)
	parameter (NStage=10,NSpan=11,NStream=11)
	common/com01/cylindrical_r(0:2*NStage*NStream-1,0:NSpan-1)
	common/com02/cylindrical_z(0:2*NStage*NStream-1,0:NSpan-1)
!	����������㾶�����������
	common/com04/velocity_axial(0:2*NStage*NStream-1,0:NSpan-1)
	common/com05/velocity_radial(0:2*NStage*NStream-1,0:NSpan-1)
	common/com06/velocity_circumferential(0:2*NStage*NStream-1,
     &	0:NSpan-1)
!	����������ϵ���������ٶ�
	common/com07/density_vapor(0:2*NStage*NStream-1,0:NSpan-1)
	common/com08/density_liquid(0:2*NStage*NStream-1,0:NSpan-1)
!	Һ��/�����ܶ�
	common/com18/volume_fraction(0:2*NStage*NStream-1,0:NSpan-1)
!	Һ���������
	common/com30/v_d_a(0:2*NStage*NStream-1,0:NSpan-1)
	common/com31/v_d_r(0:2*NStage*NStream-1,0:NSpan-1)
	common/com32/v_d_c(0:2*NStage*NStream-1,0:NSpan-1)
!	һ��Һ��������ϵ���������ٶ�
	double precision slr(NStage,0:NSpan-1),slz(NStage,0:NSpan-1)
	double precision str(NStage,0:NSpan-1),stz(NStage,0:NSpan-1)
	double precision rlr(NStage,0:NSpan-1),rlz(NStage,0:NSpan-1)
	double precision rtr(NStage,0:NSpan-1),rtz(NStage,0:NSpan-1)
!	���������Ͼ�Ҷ�Ͷ�Ҷǰ��Ե����,������ֵ����
	double precision s_l_velocity(0:NStream-1)
	double precision r_l_velocity(0:NStream-1)
!	����Ҷ/��Ҷͨ����Һ���ٶ�
	double precision s_v_velocity(0:NStream-1)
	double precision r_v_velocity(0:NStream-1)
!	����Ҷ/��Ҷͨ���������ٶ�
	double precision s_l_density(0:NStream-1)
	double precision r_l_density(0:NStream-1)
!	����Ҷ/��Ҷͨ����Һ���ܶ�
	double precision s_v_density(0:NStream-1)
	double precision r_v_density(0:NStream-1)
!	����Ҷ/��Ҷͨ���������ܶ�
	double precision s_fraction(0:NStream-1)
	double precision r_fraction(0:NStream-1)
!	����Ҷ/��Ҷͨ��Һ���������
	double precision s_weizhi(0:NStream-1)
	double precision r_weizhi(0:NStream-1)
!	����Ҷ/��Ҷͨ����λ��
	double precision area_stator(NStage,0:NSpan-1)
	double precision area_rotor(NStage,0:NSpan-1)
!	������Ҷ/��ҶβԵ���������
	double precision l_velocity_stator(NStage,0:NSpan-1)
	double precision l_velocity_rotor(NStage,0:NSpan-1)
!	������Ҷ/��ҶβԵҺ�������ٶ�
	double precision v_velocity_stator(NStage,0:NSpan-1)
	double precision v_velocity_rotor(NStage,0:NSpan-1)
!	������Ҷ/��ҶβԵ���������ٶ�
	double precision l_density_stator(NStage,0:NSpan-1)
	double precision l_density_rotor(NStage,0:NSpan-1)
!	������Ҷ/��ҶβԵ��Һ���ܶ�
	double precision v_density_stator(NStage,0:NSpan-1)
	double precision v_density_rotor(NStage,0:NSpan-1)
!	������Ҷ/��ҶβԵ�������ܶ�
	double precision fraction_stator(NStage,0:NSpan-1)
	double precision fraction_rotor(NStage,0:NSpan-1)
!	������Ҷ/��ҶβԵ��Һ���������
	double precision liquid_s_massflow(NStage,0:NSpan-1)
	double precision liquid_r_massflow(NStage,0:NSpan-1)
!	������Ҷ/��Ҷ��Ҷ��Һ������
	double precision vapor_s_massflow(NStage,0:NSpan-1)
	double precision vapor_r_massflow(NStage,0:NSpan-1)
!	������Ҷ/��Ҷ��Ҷ����������
!!!!!!
	call stator_rotor_a_r(slr,slz,str,stz,rlr,rlz,rtr,rtz)
	do k=1,NStage
		area_stator(k,0)=((str(k,1)+str(k,0))*0.5)**2-str(k,0)**2
		area_stator(k,NSpan-1)=str(k,NSpan-1)**2-
     &		((str(k,NSpan-1)+str(k,NSpan-2))*0.5)**2
		area_rotor(k,0)=((rtr(k,1)+rtr(k,0))*0.5)**2-rtr(k,0)**2
		area_rotor(k,NSpan-1)=rtr(k,NSpan-1)**2-
     &		((rtr(k,NSpan-1)+rtr(k,NSpan-2))*0.5)**2
		do j=1,NSpan-2
			area_stator(k,j)=((str(k,j+1)+str(k,j))*0.5)**2-
     &			((str(k,j)+str(k,j-1))*0.5)**2
			area_rotor(k,j)=((rtr(k,j+1)+rtr(k,j))*0.5)**2-
     &			((rtr(k,j)+rtr(k,j-1))*0.5)**2
		enddo
		do j=0,NSpan-1
		area_stator(k,j)=area_stator(k,j)*3.1415926
		area_rotor(k,j)=area_rotor(k,j)*3.1415926
		enddo
	enddo
	do k=1,NStage
		do j=0,NSpan-1
			do i=0,NStream-1
				s_weizhi(i)=cylindrical_z(2*(k-1)*NStream+i,j)
				s_l_velocity(i)=v_d_a(2*(k-1)*NStream+i,j)
				s_v_velocity(i)=velocity_axial(2*(k-1)*NStream+i,j)
				s_l_density(i)=density_liquid(2*(k-1)*NStream+i,j)
				s_v_density(i)=density_vapor(2*(k-1)*NStream+i,j)
				s_fraction(i)=volume_fraction(2*(k-1)*NStream+i,j)

				r_weizhi(i)=cylindrical_z((2*k-1)*NStream+i,j)
				r_l_velocity(i)=v_d_a((2*k-1)*NStream+i,j)
				r_v_velocity(i)=velocity_axial((2*k-1)*NStream+i,j)
				r_l_density(i)=density_liquid((2*k-1)*NStream+i,j)
				r_v_density(i)=density_vapor((2*k-1)*NStream+i,j)
				r_fraction(i)=volume_fraction((2*k-1)*NStream+i,j)
			enddo
			call SPLINE(s_weizhi,s_l_velocity,NStream-1,
     &			stz(k,j),l_velocity_stator(k,j),temp)
			call SPLINE(s_weizhi,s_v_velocity,NStream-1,
     &			stz(k,j),v_velocity_stator(k,j),temp)
			call SPLINE(s_weizhi,s_l_density,NStream-1,
     &			stz(k,j),l_density_stator(k,j),temp)
			call SPLINE(s_weizhi,s_v_density,NStream-1,
     &			stz(k,j),v_density_stator(k,j),temp)
			call SPLINE(s_weizhi,s_fraction,NStream-1,
     &			stz(k,j),fraction_stator(k,j),temp)

    			call SPLINE(r_weizhi,r_l_velocity,NStream-1,
     &			rtz(k,j),l_velocity_rotor(k,j),temp) 
     			call SPLINE(r_weizhi,r_v_velocity,NStream-1,
     &			rtz(k,j),v_velocity_rotor(k,j),temp)  	
    			call SPLINE(r_weizhi,r_l_density,NStream-1,
     &			rtz(k,j),l_density_rotor(k,j),temp)
    			call SPLINE(r_weizhi,r_v_density,NStream-1,
     &			rtz(k,j),v_density_rotor(k,j),temp)
    			call SPLINE(r_weizhi,r_fraction,NStream-1,
     &			rtz(k,j),fraction_rotor(k,j),temp)
		enddo
	enddo
	do k=1,NStage
		do j=0,NSpan-1
	liquid_s_massflow(k,j)=area_stator(k,j)*l_velocity_stator(k,j)*
     &	l_density_stator(k,j)*fraction_stator(k,j)
	liquid_r_massflow(k,j)=area_rotor(k,j)*l_velocity_rotor(k,j)*
     &	l_density_rotor(k,j)*fraction_rotor(k,j)
	vapor_s_massflow(k,j)=area_stator(k,j)*v_velocity_stator(k,j)*
     &	v_density_stator(k,j)*(1.0-fraction_stator(k,j))
	vapor_r_massflow(k,j)=area_rotor(k,j)*v_velocity_rotor(k,j)*
     &	v_density_rotor(k,j)*(1.0-fraction_rotor(k,j))
		enddo
	enddo
c	�����������������Щ����Ҫ����!
!!!!!!
	END
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	���������ٶ�
	SUBROUTINE fluctuating_velocity(f_v)
	implicit double precision (a-h,o-z)
	parameter (NStage=10,NSpan=11,NStream=11)
	common/com21/turbulence_kinetic_energy(0:2*NStage*NStream-1,
     &	0:NSpan-1)
	double precision lbound,rbound,len,my_random,kys
!	��������������
	double precision f_v(0:2*NStage*NStream-1,0:NSpan-1)
!	�����ٶ�
	lbound=-1.0
	rbound=1.0
	len=rbound-lbound
	call random_seed()
!	�����������
	do i=0,2*NStage*NStream-1
		do j=0,NSpan-1
			call random_number(my_random)
!	����0.0��1.0֮��������
			kys=len*my_random+lbound
!	����������䵽-1.0��1.0֮��
			f_v(i,j)=kys*sqrt(2/3.0*turbulence_kinetic_energy(i,j))
c			f_v(i,j)=0.0
!	���������ٶȣ���Ϊ0������Ϊû��������ɢ����
		enddo
	enddo
	END
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	�������������Ὺʼ������λ��
!	��ˮ��ֱ��С��1D-10,����Ϊû�з�������
	SUBROUTINE condensate_start(bd)
	implicit double precision (a-h,o-z)
	parameter (NStage=10,NSpan=11,NStream=11)
	common/com09/particle_diameter(0:2*NStage*NStream-1,0:NSpan-1)
	integer bd(0:NSpan-1)
	do i=0,2*NStage*NStream-2
	if (particle_diameter(i,0)<1D-10) bd(0)=i+1
	if (particle_diameter(i,NSpan-1)<1D-10) bd(NSpan-1)=i+1
		do j=1,NSpan-2
			if (particle_diameter(i,j+1)<1D-10.or.
     &			particle_diameter(i,j)<1D-10.or.
     &			particle_diameter(i,j-1)<1D-10.or.
     &			particle_diameter(i+1,j)<1D-10) then
				bd(j)=i+1
			end if
		enddo
	enddo
	END
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	��������ѧ��ʧ����
!	�ο���Fortran77 �ṹ��������ơ���̷��ǿ�ȣ�P176
	SUBROUTINE thermodynamicloss()
	implicit double precision (a-h,o-z)
	parameter (NStage=10,NSpan=11,NStream=11)
	common/com16/axial_distance(0:2*NStage*NStream-1)
	double precision at,bt
	integer ns
	double precision loss(1:2*NStage)
	double precision thermodynamic_loss(NStage)
	external thermo_loss
	do i=1,2*NStage
	at=axial_distance(NStream*(i-1))
	bt=axial_distance(NStream*i-1)
	ns=1
	call fsim(at,bt,ns,thermo_loss,loss(i))
	enddo
	do i=1,NStage
		thermodynamic_loss(i)=loss(2*i-1)+loss(2*i)
		if (thermodynamic_loss(i)<=1.0) thermodynamic_loss(i)=0.0
	enddo
	print *,'����ѧ��ʧ:'
	do i=1,NStage
		print *,i,thermodynamic_loss(i)
	enddo
!	����ѧ��ʧ�洢��thermodynamic_loss�У����������洢
	END
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	����Һ�ζ��μд���ʧ����
	SUBROUTINE reentainmentloss()
	implicit double precision (a-h,o-z)
	parameter (NStage=10,NSpan=11,NStream=11)
	common/com02/cylindrical_z(0:2*NStage*NStream-1,0:NSpan-1)
	double precision ad,bd
!	ad:����������,bd:����������
	integer L
!	�������������ѡ��
	double precision re_entainment_a_loss(1:2*NStage)
	double precision re_entainment_c_loss(1:2*NStage)
	double precision re_entainment_r_loss(1:2*NStage)
	double precision loss(1:2*NStage)
	double precision re_entainment_loss(NStage)
!	���μд���ʧ
	external D_L_BOUNDARY,DRAG_LOSS
	do i=1,2*NStage
	print *,'i=',i
	ad=(cylindrical_z(NStream*(i-1),0)+
     &	cylindrical_z(NStream*(i-1),NSpan-1))*0.5
	bd=(cylindrical_z(NStream*i-1,0)+
     &	cylindrical_z(NStream*i-1,NSpan-1))*0.5
	L=1
	call FSIM2(ad,bd,D_L_BOUNDARY,DRAG_LOSS,L,re_entainment_a_loss(i))
	L=2
	call FSIM2(ad,bd,D_L_BOUNDARY,DRAG_LOSS,L,re_entainment_c_loss(i))
	L=3
	call FSIM2(ad,bd,D_L_BOUNDARY,DRAG_LOSS,L,re_entainment_r_loss(i))
	enddo
	do i=1,2*NStage
		loss(i)=re_entainment_a_loss(i)+re_entainment_c_loss(i)+
     &		re_entainment_r_loss(i)
	enddo
	do i=1,NStage
		re_entainment_loss(i)=loss(2*i-1)+loss(2*i)
		if (re_entainment_loss(i)<=1.0) re_entainment_loss(i)=0.0
	enddo
	print *,'���μд���ʧ:'
	do i=1,NStage
		print *,i,re_entainment_loss(i)
	enddo
!	Һ�ζ��μд���ʧ�洢��re_entainment_loss�У����������洢
	END
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	�����ƶ���ʧ����
	SUBROUTINE impactloss()
	parameter (NStage=10,NSpan=11,NStream=11)
	implicit double precision (a-h,o-z)
	double precision slr(NStage,0:NSpan-1),slz(NStage,0:NSpan-1)
	double precision str(NStage,0:NSpan-1),stz(NStage,0:NSpan-1)
	double precision rlr(NStage,0:NSpan-1),rlz(NStage,0:NSpan-1)
	double precision rtr(NStage,0:NSpan-1),rtz(NStage,0:NSpan-1)
!	��Ҷ�Ͷ�Ҷǰ��Ե�������������
	double precision r1,r2
	double precision impact_loss(NStage)
	external DEOPLET_IMPACT
	call stator_rotor_a_r(slr,slz,str,stz,rlr,rlz,rtr,rtz)	
	do k=1,NStage
	print *,'k=',k
		r1=rlr(k,0)
		r2=rlr(k,NSpan-1)
		call fsim(r1,r2,k,DEOPLET_IMPACT,impact_loss(k))
		if (impact_loss(k)<=1.0) impact_loss(k)=0.0
	enddo
	print *,'�ƶ���ʧ:'
	do i=1,NStage
		print *,i,impact_loss(i)
	enddo
!	����Һ���ƶ���ʧ�洢��impact_loss�У����������洢
	END
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	������ˮ��ʧ����
	SUBROUTINE collectedwaterloss()
	implicit double precision (a-h,o-z)
	parameter (NStage=10,NSpan=11,NStream=11)
	common/com80/turbine_power(NStage)
!	��������
	double precision dp_rt_s(NStage,0:NSpan-1)
	double precision dp_rt_r(NStage,0:NSpan-1)
!	һ��Һ���ھ�Ҷ�Ͷ�Ҷ�ĳ�����
	double precision liquid_s_massflow(NStage,0:NSpan-1)
	double precision liquid_r_massflow(NStage,0:NSpan-1)
!	һ��Һ���ھ�Ҷ�Ͷ�Ҷ������
	double precision vapor_s_massflow(NStage,0:NSpan-1)
	double precision vapor_r_massflow(NStage,0:NSpan-1)
!	�����ھ�Ҷ�Ͷ�Ҷ������
	double precision deposition_s_massflow(NStage,0:NSpan-1)
	double precision deposition_r_massflow(NStage,0:NSpan-1)
!	һ��Һ���ڸ�����Ҷ�ߵĳ�����
	double precision liquid_d_massflow(NStage)
	double precision massflow(NStage)
!	һ�γ���Һ�κ�ʪ�����������ڸ����ķֲ�
	double precision temp
	double precision collectedwater_loss(NStage)
	call DEPOSITION_RATE(dp_rt_s,dp_rt_r)
!	���������
	call FLOWRATE(liquid_s_massflow,liquid_r_massflow,
     &	vapor_s_massflow,vapor_r_massflow)
!	����Һ��������ڸ�����Ҷ�ߵ�����
	do i=1,NStage
		do j=0,NSpan-1
	deposition_s_massflow(i,j)=dp_rt_s(i,j)*liquid_s_massflow(i,j)
	deposition_r_massflow(i,j)=dp_rt_r(i,j)*liquid_r_massflow(i,j)
		enddo
	enddo
!	����һ��Һ���ڸ�����Ҷ�Ͷ�Ҷ��Ҷ�ߵĳ�����
	do i=1,NStage
		liquid_d_massflow(i)=0.0
		massflow(i)=0.0
		do j=0,NSpan-1
			liquid_d_massflow(i)=liquid_d_massflow(i)+
     &			deposition_s_massflow(i,j)+deposition_r_massflow(i,j)
			massflow(i)=massflow(i)+liquid_r_massflow(i,j)+
     &			vapor_r_massflow(i,j)
		enddo
	enddo
!	����һ�γ���Һ�κ������ڸ���������
	do i=1,NStage
		collectedwater_loss(i)=liquid_d_massflow(i)/massflow(i)*
     &		turbine_power(i)
		if (collectedwater_loss(i)<=1.0) collectedwater_loss(i)=0.0
	enddo
	print *,'��ˮ��ʧ:'
	do i=1,NStage
		print *,i,collectedwater_loss(i)
	enddo
!	��ˮ��ʧ�洢��collectedwater_loss�У����������洢
	END
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	����������ʧ����
	SUBROUTINE centrifugingloss()
	implicit double precision (a-h,o-z)
	parameter (NStage=10,NSpan=11,NStream=11)
	common/com03/omiga(0:2*NStage*NStream-1,0:NSpan-1)
	double precision dp_rt_s(NStage,0:NSpan-1)
	double precision dp_rt_r(NStage,0:NSpan-1)
!	һ��Һ���ھ�Ҷ�Ͷ�Ҷ�ĳ�����
	double precision liquid_s_massflow(NStage,0:NSpan-1)
	double precision liquid_r_massflow(NStage,0:NSpan-1)
!	һ��Һ���ھ�Ҷ�Ͷ�Ҷ������
	double precision vapor_s_massflow(NStage,0:NSpan-1)
	double precision vapor_r_massflow(NStage,0:NSpan-1)
!	�����ھ�Ҷ�Ͷ�Ҷ������
	double precision deposition_r_massflow(NStage,0:NSpan-1)
!	һ��Һ���ڶ�Ҷ��Ҷ�ߵĳ�����
	double precision massflow_r_deposition(NStage)
!	һ��Һ���ڸ�����Ҷ������
	double precision slr(NStage,0:NSpan-1),slz(NStage,0:NSpan-1)
	double precision str(NStage,0:NSpan-1),stz(NStage,0:NSpan-1)
	double precision rlr(NStage,0:NSpan-1),rlz(NStage,0:NSpan-1)
	double precision rtr(NStage,0:NSpan-1),rtz(NStage,0:NSpan-1)
!	��Ҷ�Ͷ�Ҷǰ��Ե�������������
	double precision centrifuging_loss(NStage)
	call DEPOSITION_RATE(dp_rt_s,dp_rt_r)
!	���������
	call FLOWRATE(liquid_s_massflow,liquid_r_massflow,
     &	vapor_s_massflow,vapor_r_massflow)
!	����Һ��������ڸ�����Ҷ�ߵ�����
	call stator_rotor_a_r(slr,slz,str,stz,rlr,rlz,rtr,rtz)
!	���������Ҷ�Ͷ�Ҷ��ǰ��Ե����
	do i=1,NStage
		do j=0,NSpan-1
	deposition_r_massflow(i,j)=dp_rt_r(i,j)*liquid_r_massflow(i,j)
		enddo
	enddo
!	����һ��Һ���ڸ�����Ҷ�Ͷ�Ҷ��Ҷ�ߵĳ�����
	do i=1,NStage
		massflow_r_deposition(i)=0.0
		do j=0,NSpan-1
			massflow_r_deposition(i)=massflow_r_deposition(i)+
     &			deposition_r_massflow(i,j)
		enddo
	enddo
!	����һ��Һ���ڸ�����Ҷ�Ͷ�Ҷ�ܳ�����
	do i=1,NStage
		centrifuging_loss(i)=massflow_r_deposition(i)*
     &		(omiga((2*i-1)*NStream,0)**2)*(((rlr(i,NSpan-1)+
     &	rtr(i,NSpan-1))*0.5)**2-((rlr(i,0)+rtr(i,0))*0.5)**2)
		if (centrifuging_loss(i)<=1.0) centrifuging_loss(i)=0.0
	enddo

	print *,'������ʧ:'
	do i=1,NStage
		print *,i,centrifuging_loss(i)
	enddo
!	������ʧ�洢��centrifuging_loss�У����������洢
	END
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
