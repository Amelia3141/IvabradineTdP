import os
import sys
import logging
import importlib

# Configure logging
logging.basicConfig(
    level=logging.DEBUG,  # Set to DEBUG for more detailed logs
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.StreamHandler(sys.stdout)
    ]
)
logger = logging.getLogger(__name__)

def check_imports():
    """Check and log import status of required modules"""
    required_modules = {
        'flask': None,
        'torch': None,
        'transformers': None,
        'sentence_transformers': None,
        'Bio': None
    }
    
    for module_name in required_modules:
        try:
            module = importlib.import_module(module_name)
            required_modules[module_name] = getattr(module, '__version__', 'unknown')
            logger.info(f"Successfully imported {module_name} (version: {required_modules[module_name]})")
        except ImportError as e:
            logger.error(f"Failed to import {module_name}: {str(e)}")
            required_modules[module_name] = f"error: {str(e)}"
    
    return required_modules

def initialize_app():
    """Initialize the application and verify the environment"""
    try:
        # Log environment information
        logger.info("Initializing application...")
        logger.info(f"Python version: {sys.version}")
        logger.info(f"Current working directory: {os.getcwd()}")
        
        # Check environment variables
        port = os.environ.get('PORT')
        logger.info(f"PORT environment variable: {port}")
        
        # Check imports first
        module_status = check_imports()
        
        # Only import app if critical modules are available
        if all(v != None for v in module_status.values()):
            from api_server import app
            logger.info("Successfully imported api_server")
            return app
        else:
            raise ImportError("Failed to import required modules")
            
    except Exception as e:
        logger.error(f"Failed to initialize application: {str(e)}")
        raise

app = None

if __name__ == "__main__":
    if "--test-only" in sys.argv:
        try:
            app = initialize_app()
            logger.info("Application initialized successfully in test mode")
            sys.exit(0)
        except Exception as e:
            logger.error(f"Test initialization failed: {str(e)}")
            sys.exit(1)
    else:
        app = initialize_app()
        port = int(os.environ.get('PORT', 10000))
        logger.info(f"Starting server on port {port}")
        app.run(host='0.0.0.0', port=port)
else:
    # For gunicorn
    logger.info("Initializing application for gunicorn")
    app = initialize_app()
