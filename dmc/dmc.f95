module dmc
    use normrand ! pro generovani nahodnych cisel z normalniho rozdeleni

    implicit none

    contains
    
    subroutine simulace(p_bloku, p_kroku, N0, N_max, N_min, alfa, dt, x0, Er, potencial)
        implicit none
        integer, intent(in) :: p_bloku, p_kroku, N0, N_max, N_min
        real, intent(in) :: alfa, dt, x0(:)
        real, intent(out) :: Er
        integer :: h, m(N_max), pra, p_0, p_replik, j, i_blok, i_krok, dp, I_0(N_max), ub ! promenne vysvetleny nize
        real :: rep(N_max, size(x0)), pe, pot, potencial
        
        m = [(1, j=1,N_max)] ! 'spawn number', udava zda-li replika zije a zda-li se ma zkopirovat
        write(*,"(1x,'Blok',4x,'N',5x,'<H>',6x,'<V>')")
        dp = size(x0) ! prostorova dimenze problemu
        rep = reshape([(x0, j=1,N_max)],[N_max,dp]) ! repliky jsou ulozeny v matici (prvni index - cislo repliky, druhy - prostorava dimenze)
        p_replik = N0 ! pocet zivych replik v danny moment
        h = N0 ! protoze ve fortranu nejde snadno za chodu menit rozmer vektoru a matic, jsou vsechny delsi nez je treba a bere se z nich jen zacetek (od 1 do h)
        Er = potencial(x0) ! inicializace referencni energie
        do i_blok = 1,p_bloku ! simulace probiha v blocich, ve kterych se sbiraji ruzne statistiky (zatim jen potencial a Er)
            pot = 0
            do i_krok = 1,p_kroku
                call difuze(rep, m, I_0, pra, p_0, pe, h, dp, Er, dt, ub, potencial) ! upravi polohu replik (a)
                call vetveni(rep, m, I_0, pra, p_0, h, Er, alfa, p_replik, ub) ! patricne pozkopiruje repliky
                pot = pot + pe ! zprumeruje potencialni enrgii v replik v bloku
            end do
            pot = pot / p_kroku
            write(*,"(I3,3x,I7,f9.5,f9.5)") i_blok, p_replik, Er, pot ! vypise udaje o probehnuvsim kroku
            ! call renormalizace()
        end do
        call uloz(rep, m, h) ! zapise polohy zivych replik do souboru
    end subroutine simulace
    
    ! =================================================
    
    subroutine difuze(rep, m, I_0, pra, p_0, pe, h, dp, Er, dt, ub, potencial)
        implicit none
        integer, intent(inout) :: m(:), I_0(:)
        integer, intent(in) :: h, dp
        real, intent(inout) :: rep(:,:)
        real, intent(out) :: pe
        real, intent(in) :: Er, dt
        integer, intent(out) :: pra, p_0, ub
        integer :: j, k
        real :: potencial, u, pot_e = 0
        
        pe = 0
        pra = 0
        k = 1
        ub = 0
        do j = 1,h ! v pouzivanem rozsahu
            if (m(j) /= 0) then ! pro vesechny zive repliky
                call random_number(u)
                rep(j,:) = rep(j,:) + sqrt(dt) * normrnd(dp) ! posune repliku
                pot_e = potencial([rep(j,:)]) ! vypocte potencial danne repliky
                pe = pe + pot_e
                m(j) = min(3, floor(u + exp((Er - pot_e)*dt))) ! vypocte odpovidajici 'spawn number'
                if (m(j) /= 0) then ! pokud neni nulovy (tzn. replika bude zit dal) 
                    pra = pra + m(j) - 1 ! 'pra' je pocet vytvorenych kopii
                else ! pokud je nulovy (tzn. replika 'zemre')
                    I_0(k) = j ! zapise si polohu repliky do seznamu volnych mist, ktere poto zaplni novymi replikami
                    k = k + 1
                    ub = ub + 1 ! 'ub' je pocet replik, ktere behem iterace 'zemrou'
                end if
            else ! pokud byla replika uz mrtva
                I_0(k) = j ! prida jeji misto na seznam cekatelu k zaplneni
                k = k + 1
            end if
        end do
        p_0 = k - 1 ! 'p_0' je pocet mrtvych replik v pouzivanem rozsahu (tzn. od 1 do h)
        pe = pe / (h - p_0) ! vypocte prumerny potencial zivych replik
    end subroutine difuze
    
    ! =================================================
    
    subroutine vetveni(rep, m, I_0, pra, p_0, h, Er, alfa, p_replik, ub)
        implicit none
        integer, intent(in) :: pra, I_0(:), p_0, ub
        real, intent(inout) :: rep(:,:), Er
        real, intent(in) :: alfa
        integer, intent(inout) :: h, p_replik, m(:)
        integer :: nove(pra), h_nove, j, k, i
    
        if (pra > p_0) then ! pokud je novych replik vice nez prazdnych mist
            nove = [I_0(1:p_0),(i, i=h+1,h+pra-p_0)] ! 'nove' je seznam mist, kam budou umisteny nove repliky
            h_nove = h + pra - p_0 ! natahne pouzivanny rozsah
        else ! pokud prazdna mista staci
            nove = I_0(1:pra)
            h_nove = h ! necha rozsah tak, jak je
        end if

        k = 1
        do j = 1,h
            if (m(j) == 2) then ! pro repliky s m = 2
                rep(nove(k),:) = rep(j,:) ! vytvori jednu kopii
                m(nove(k)) = 1 ! nastavi 'm' nove repliky na 1, jinak by byla v subrutine 'difuze' povazovana za mrtvou
                k = k + 1
            elseif (m(j) == 3) then ! pro repliky s m = 3
                rep(nove([k,k+1]),:) = rep([j,j],:) ! vytvori 2 kopie
                m(nove([k,k+1])) = [1,1]
                k = k + 2
            endif
        end do
        
        ! jestli h_nove neni v <N_min, N_max>, renormalizuje
        !if (h_nove > N_max) then
            !renormalizace a konec bloku
        !elseif (h_nove < N_min) then
            !renormalizace a konec bloku
        !else
            h = h_nove ! upravi pouzivanny rozsah
            Er = Er + alfa * (1 - real(p_replik + pra - ub) / p_replik) ! upravi referencni energii
            p_replik = p_replik + pra - ub ! upravi pocet replik
        !endif
    
    end subroutine vetveni

    subroutine uloz(rep, m, h) ! zapise konecnou populaci zivych replik do souboru
        implicit none
        integer :: j
        integer, intent(in) :: h, m(:)
        real, intent(in) :: rep(:,:)
        open(1,file='populace_dmc')
        do j=1,h
            if (m(j) /= 0) then
                write(1,*) rep(j,:)
            end if
        end do
        close (1)
    end subroutine uloz

end module dmc
