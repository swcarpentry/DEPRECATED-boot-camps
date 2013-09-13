from morse import MorseTranslator

class TestMorseTranslator:

        def test(self):
            translator = MorseTranslator()
            assert "... --- ..." == translator.encode("SOS")
            assert "sos" == translator.decode("... --- ...")
            print "OK"

    	if __name__ == "__main__":    

        	test_translator = TestMorseTranslator()
        	test_translator.test()