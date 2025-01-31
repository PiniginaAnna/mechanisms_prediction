pseudoatoms_dict = {'Nu-:': ['[F-;X0]',
                             '[Cl-;X0]',
                             '[Br-;X0]',
                             '[I-;X0]',
                             '[N-;X2]',
                             '[N-;X1]',
                             '[O-;X1]',
                             '[S-;X1]',
                             '[C-;X1,X3]',
                             # '[C-;X1](#[N+0:100])',
                             # '[O-;X1](-[N+0;X2:100]=[O+0:101])'
                             '[#1-]'],
                    'Nu+0': ['[O+0;X1]',
                             '[S+0;X1]',
                             '[N+0;X3]',
                             '[N+0;X2]',
                             '[n+0;X2]'
                             '[P+0;X3]',
                             '[O+0;X2]',
                             # '[N+0;X2](=[O+0:100])(-[O-:101])',
                             '[S+0;X2]'],
                    'E+:': ['[C+;X3,X2,X1]', '[#1+]', '[N+;X2](=[O+0:100])(=[O+0:101])'],
                    'E+0': ['[F+0;X1]',
                            '[Cl+0;X1]',
                            '[Br+0;X1]',
                            '[I+0;X1]',
                            '[S+0;X3](=[O+0:100])(=[O+0:101])(=[O+0:102])',
                            '[S+0;X3](=[O+;H1:100])(=[O+0:101])(=[O+0:102])'],
                    'L_nuc+0': ['[F+0;X1]',
                                '[Cl+0;X1]',
                                '[Br+0;X1]',
                                '[I+0;X1]',
                                '[N+0;X3]([#1+0:100])([#1+0:101])',
                                '[N+0](-[C:100]=[O:101])',
                                '[N+0;X3](-[C:100](=[O:101])-[O:102])',
                                '[N+0](-[C:100](=[O:101])-[N:102])',
                                '[C+0](=[O:100])(-[O+0:101]-[C:102])',
                                '[C+0;X4]([F+0,Cl+0,Br+0,I+0;X1:100])([F+0,Cl+0,Br+0,I+0;X1:101])([F+0,Cl+0,Br+0,I+0;X1:102])',
                                '[O+0;X2]([P+0:100]([F+0,Cl+0,Br+0,I+0;X1:101])([F+0,Cl+0,Br+0,I+0;X1:102]))',
                                '[O+0;X2]([P+0:100](=[O+0:101])([F+0,Cl+0,Br+0,I+0;X1:102])[F+0,Cl+0,Br+0,I+0;X1:103])',
                                '[O+0;X2]([P+0:100](=[O+0:101])(=[O+0:102]))',
                                '[O+0;X2]([P+:100](=[O+0:101])(-[O+0;H1:102]))',
                                '[O+0;X2]([N:100]1[N:101]=[N:102][C:103]2=[C:104][C:105]=[C:106][C:107]=[C:108]12)',
                                '[O+0;X2]([N:100]1[N:101]=[N:102][C:103]2=[C:104][C:105]=[C:106][N:107]=[C:108]12)',
                                '[O+0;X2]([C+0:100](=[N+:101]([#1+0:102])[C+0:103])[N+0:104]([#1+0:105])[C+0:106])',
                                '[O+0;X2]([C+0:100](=[O+0:101])[C+0:102])',
                                '[O+0;X2](-[S:100](=[O:101])[O:102][C:103])',
                                '[O+0;X2](-[S:100](=[O:101])(=[O:102])[C:103])',
                                '[O+0;X2](-[C:100](=[O:101])-[N:102])',
                                '[O+0;X2]([#1+0:100])',
                                '[O+0;X2](-[S:100](=[O:101])(=[O:102])[C:103]1=[C:104][C:105]=[C:106]([C:107])[C:108]=[C:109]1)',
                                '[O+0;X2]([S:100](=[O:101])(=[O:102])[C:103]1=[C:104][C:105]=[C:106]([Br:107])[C:108]=[C:109]1)',
                                '[O+0;X2]([S:100](=[O:101])(=[O:102])[C:103]1=[C:104][C:105]=[C:106]([N+:107](=[O:108])[O-:109])[C:110]=[C:111]1)',
                                '[O+0;X2](-[S:100](=[O:101])([C:102])=[O:103])',
                                '[O+0;X2]([S:100](=[O:101])(=[O:102])[C:103]([F:104])([F:105])[F:106])',
                                '[O+0;X2]([C;X4,X3:100])',
                                '[O+0;X2](-[C:100](=[O:101])-[O:102]-[C:103])',
                                '[O+0;X2](-[P+0;X4:100])',
                                '[S+0;X2](-[C,#1;+0:100])',
                                '[S+0](=[O:100])(=[O:101])([C:102]1=[C:103][C:104]=[C:105]([C:106])[C:107]=[C:108]1)',
                                '[S+0;X2](-[C+0:100]#[N+0:101])',
                                '[S+0;X4](=[O+0:100])(=[O+0:101])([C+0:102])',
                                '[O+0;X2](-[S:100](=[O:101])[Cl:102])'],
                    'L+:': ['[O+;X3]([#1+0:100])([#1+0:101])',
                            '[N+;X4]',
                            '[N+;X2](#[N+0:100])',
                            '[O+;X2](=[P+0:100]([C+0:101])([C+0:102])([C+0:103]))',
                            '[N+;X3](=[O+0:100])(-[#1+0:101])',
                            '[N+;X3](=[O+0:100])(-[O-:101])',
                            '[S+](=[O+0:100])([C+0:101])([C+0:102])',
                            '[S+]([C+0:100])([C+0:101])',
                            '[O+;X3]([C;X4,X3:100])',
                            '[O+;X3]([P+0:100]([F+0,Cl+0,Br+0,I+0;X1:101])([F+0,Cl+0,Br+0,I+0;X1:102]))',
                            '[P+;H3;X4]'],
                    'L_el+0': ['[#1+0]'],
                    'Bs-:': ['[O-;X1]', '[C-;X3,X1]', '[#1-]', '[N-;X2]', '[S-;X1]'],
                    'Bs+0': ['[O+0;X2]([#1:100])([#1:101])',
                             '[N+0;X2,X3]',
                             '[O+0;X2]',
                             '[S+0;X2]',
                             '[C0+;X4]([C:100][C:101][C:102])([Li:103])'],
                    'Ac+0': ['[O+0;H1,H2;X2]',
                             '[N+0;H1,H2,H3;X3]',
                             '[F+0;H1]',
                             '[Cl+0;H1]',
                             '[Br+0;H1]',
                             '[I+0;H1]'],
                    'Ac+:': ['[O+;H3;X3]([#1:100])([#1:101])',
                             '[N+;H1,H2,H3,H4;X3,X4]',
                             '[O+;H1,H2;X3]'],
                    'LA+0': ['[Al+0;X3]', '[Fe+0;X3]'],
                    'LA-:': ['[Al-;X4]', '[Fe-;X4]'],
                    'LB+0': ['[F+0;X1]',
                             '[Cl+0;X1]',
                             '[Br+0;X1]',
                             '[I+0;X1]',
                             '[O+0;X1]',
                             '[O+0;X2]'],
                    'LB+:': ['[F+;X2]', '[Cl+;X2]', '[Br+;X2]', '[I+;X2]', '[O+;X2]'],
                    'EWG+0': ['[N+;X2](#[N+0:100])',
                              '[C+0](=[O+0:100])([C:101])',
                              '[C+0;X3](=[O+;H1:100])',
                              '[C+0;X2](#[N+0:100])',
                              '[C+0;X3](-[O+0:100])(=[O+0:101])',
                              '[S+](=[O+0:100])([C+0:101])([C+0:102])',
                              '[S+](-[O+0:100])',
                              '[S+0](-[O+0:100])',
                              '[O+0](-[S+:100])',
                              '[P+]([C+0:100])([C+0:101])([C+0:102])',
                              '[C+0](=[O+0:100])([F+0,Cl+0,Br+0,I+0;X1:101])',
                              '[C+0](=[O+0:100])(-[#1:101])',
                              '[C+0](=[O+0:100])(-[N;X3:101])',
                              # '[C+0](=[O+0:100])(-[O+0:101]-[C]=[O+0:102])',
                              '[C+0](=[O:100])(-[O+0:101]-[C:102])',
                              '[C+0](=[N:100]-[C,#1:101])',
                              '[C+0](=[N:100]-[O+0:101]-[#1:102])',
                              '[C+0](=[O+0:100])(-[S:101]-[C:102])',
                              '[C+0;X4]([F+0,Cl+0,Br+0,I+0;X1:100])([F+0,Cl+0,Br+0,I+0;X1:101])([F+0,Cl+0,Br+0,I+0;X1:102])',
                              '[C+0](-[N+:100](=[O+0:101])-[O-:102])',
                              '[C+0](-[S:100](=[O+0:101])(=[O+0:102])-[O+0:103])',
                              '[C+0](-[S:100](=[O+0:101])(=[O+0:102])-[C:103])',
                              '[C+0](-[S:100](=[O+0:101])-[C:102])'],
                    'EDG+0': []}
