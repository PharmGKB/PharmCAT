setup:
	bundle config set --local path 'vendor/bundle'
	bundle install

serve:
	bundle exec jekyll serve --incremental

build:
	bundle exec just-the-docs rake search:init
	bundle exec jekyll build

clean:
	rm -rf _site

next: clean build
	surge _site next.pharmcat.org

dev:
	bundle exec just-the-docs rake search:init
	bundle exec jekyll build --incremental
	cd ..; yarn localDocs

local:
	bundle exec jekyll serve --force-polling --incremental
