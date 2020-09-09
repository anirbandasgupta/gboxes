"""frame URL Configuration

The `urlpatterns` list routes URLs to views. For more information please see:
    https://docs.djangoproject.com/en/3.0/topics/http/urls/
Examples:
Function views
    1. Add an import:  from my_app import views
    2. Add a URL to urlpatterns:  path('', views.home, name='home')
Class-based views
    1. Add an import:  from other_app.views import Home
    2. Add a URL to urlpatterns:  path('', Home.as_view(), name='home')
Including another URLconf
    1. Import the include() function: from django.urls import include, path
    2. Add a URL to urlpatterns:  path('blog/', include('blog.urls'))
"""
from django.contrib import admin
from django.urls import path
from django.conf import settings
from django.conf.urls.static import static
from input01 import views

urlpatterns = [
    path('admin/', admin.site.urls),
    path('', views.welcome, name="welcome"),
    path('ourResponse', views.ourResponse, name="ourResponse"),
    path('takeInput', views.takeInput, name="takeInput"),
    path('performAlgo',views.performAlgo, name="performAlgo"),
    path('downloadfile',views.downloadfile, name="downloadfile"),
    path('downloadfilenew',views.downloadfilenew, name="downloadfilenew"),
    path('negative_control',views.negative_control, name="negative_control"),
    path('neg_control',views.neg_control, name="neg_control"),
    path('negative1',views.negative1, name="negative1"),
    path('negative2',views.negative2, name="negative2"),
    path('negative3',views.negative3, name="negative3"),
    path('downloadfile1',views.downloadfile1, name="downloadfile1"),
    path('showimage', views.showimage, name='showimage'),
    path('showimage1', views.showimage1, name='showimage1')
]

if settings.DEBUG:
    urlpatterns += static(settings.MEDIA_URL, document_root = settings.MEDIA_ROOT)