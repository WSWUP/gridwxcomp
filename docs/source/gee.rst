.. _getting_started_gee:

Getting Started with Google Earth Engine
========================================

This section will guide you through the process of setting up Google Earth Engine (GEE) to work with ``gridwxcomp``. Follow these steps to ensure you have properly set up your GEE account before proceeding with the tutorial.

Account creation
----------------

Before you can use Google Earth Engine, you need to create an account:

1. **Create a Google Earth Engine Account**:

   If you haven't already, visit the `Google Earth Engine Sign Up Page <https://earthengine.google.com/signup/>`_ to create an account.
   Approval may take a few minutes to a few days, so it's best to do this in advance.

2. **Authenticate with Google Earth Engine**:

   After installing the GEE API (should already be installed with ``gridwxcomp`` environment), authenticate with Google Earth Engine from a terminal window:

   .. code-block:: bash

      earthengine authenticate

   This command will open a browser window asking for your Google account credentials and permission to allow access. Once authenticated, you will receive a verification code that must be pasted back into the terminal to complete the setup.

Create a GEE Project
--------------------

Starting Nov. 12, 2024, `GEE requires a cloud project <https://developers.google.com/earth-engine/guides/transition_to_cloud_projects>`_ to be specified for all API access without a quota. After authenticating your account you should proceed with these steps to set up your first GEE cloud project name which only takes a few minutes, after this you can initialize GEE in your Python environment. 

**Create a Google Cloud Project for GEE**:

   - Go to the `Google Cloud Console <https://console.cloud.google.com/>`_ and sign in using your Google account.
   - Click on **"New Project"** to create a new project.
   - Give your project a name and note the **Project ID** (this ID will be used in the GEE initialization).
   - Make sure to **enable the Google Earth Engine API** for this project by navigating to the **API Library** in the Cloud Console.

Verify your GEE API account access
----------------------------------
Here are a few steps to take to ensure your installation of GEE Python API went smoothly. 

1. **Specify Project ID for Initialization**:

   - When initializing Earth Engine, provide your Project ID like this:

   .. code-block:: python

      import ee
      ee.Initialize(project='my-project-id')

   Replace **'my-project-id'** with the actual Project ID you created in Google Cloud.

   If you encounter any errors during initialization, double-check that you have authenticated properly, that the Project ID is correct, and that your environment is set up correctly. 

2. **Test Account Access**:

   You can further verify your Google Earth Engine setup by checking your user information:

   .. code-block:: python

      info = ee.data.getAssetRoots()
      print(info)

   If your account is correctly set up, this should return a list of available assets or an empty list if you haven't created any assets yet.

Once your Google Earth Engine account is set up and verified, you can proceed to the ``gridwxcomp`` :ref:`tutorial <tutorial>`.


