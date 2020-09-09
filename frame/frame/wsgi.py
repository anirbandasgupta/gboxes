"""
WSGI config for frame project.

It exposes the WSGI callable as a module-level variable named ``application``.

For more information on this file, see
https://docs.djangoproject.com/en/3.0/howto/deployment/wsgi/
"""

import os, sys

from django.core.wsgi import get_wsgi_application

sys.path.append('/var/www/html/datascience')

os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'frame.settings')

application = get_wsgi_application()
