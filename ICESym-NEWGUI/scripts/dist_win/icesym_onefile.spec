# -*- mode: python -*-

block_cipher = None

added_files = [
               (os.environ['ANACONDA3PATH']+'\\Library\\plugins\\platforms\\qwindows.dll', 'platforms')
              ]

a = Analysis(['..\\__main__.py'],
             binaries=[],
             datas=added_files,
             hiddenimports=[],
             hookspath=[],
             runtime_hooks=[],
             excludes=[],
             win_no_prefer_redirects=False,
             win_private_assemblies=False,
             cipher=block_cipher,
             noarchive=False)
pyz = PYZ(a.pure, a.zipped_data,
             cipher=block_cipher)
exe = EXE(pyz,
          a.scripts,
	     a.binaries - TOC([('mkl_avx512_mic.dll',None,None),('mkl_avx512.dll',None,None),('mkl_avx2.dll',None,None),('mkl_mc3.dll',None,None),('mkl_avx.dll',None,None),('mkl_mc.dll',None,None),('mkl_rt.dll',None,None),('mkl_vml_avx2.dll',None,None),('mkl_vml_avx.dll',None,None),('mkl_sequential.dll',None,None),('mkl_vml_avx512.dll',None,None),('mkl_vml_mc3.dll',None,None),('mkl_vml_mc.dll',None,None),('mkl_vml_mc2.dll',None,None),('mfc140u.dll',None,None),('mkl_vml_def.dll',None,None),('mkl_vml_cmpt.dll',None,None)]),
          a.zipfiles,
          a.datas,
          [],
          name='icesym',
          debug=False,
          bootloader_ignore_signals=False,
          strip=False,
          upx=True,
          runtime_tmpdir=None,
          console=False,
          icon='engine.ico' )
